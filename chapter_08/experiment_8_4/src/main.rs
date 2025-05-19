use anyhow::{Context, Result, anyhow};
use clap::Parser;
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn, error};
use noodles_vcf as vcf;
use noodles_gff as gff;
use noodles_bgzf as bgzf;
use bio::io::fasta::IndexedReader;
use polars::prelude::*;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use serde::{Serialize, Deserialize};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    sync::{Arc, Mutex},
    time::Instant,
};
// Temporarily commenting out tch imports
// use tch::{CModule, Tensor, Device};

// Define stub types to keep the code compiling
#[derive(Debug, Clone)]
struct CModule;
#[derive(Debug, Clone)]
struct Tensor;
#[derive(Debug, Clone)]
struct Device;
use thiserror::Error;

/// Errors specific to variant annotation
#[derive(Error, Debug)]
pub enum AnnotationError {
    #[error("Failed to load GFF file: {0}")]
    GffLoadError(String),
    
    #[error("Failed to load VCF file: {0}")]
    VcfLoadError(String),
    
    #[error("Failed to load frequencies from {0}: {1}")]
    FreqLoadError(String, String),
    
    #[error("Failed to load PyTorch model: {0}")]
    ModelLoadError(String),
    
    #[error("Failed to fetch sequence: {0}")]
    SequenceFetchError(String),
    
    #[error("No reference genome provided for splice prediction")]
    NoReferenceError,
    
    #[error("Failed to predict splice effect: {0}")]
    SplicePredictionError(String),
}

/// Command line arguments
#[derive(Parser, Debug, Clone)]
#[clap(
    name = "variant-annotator",
    version = "0.2.0",
    author = "Genomics Team",
    about = "Annotates genetic variants with gene, frequency, and splice effects"
)]
struct Args {
    /// Input VCF file path
    #[arg(short, long)]
    vcf: String,
    
    /// Gene annotation GFF file path
    #[arg(short, long)]
    gff: String,
    
    /// gnomAD frequency data file path (compressed)
    #[arg(short, long)]
    gnomad: String,
    
    /// Optional reference genome FASTA for splice predictions
    #[arg(short, long)]
    reference: Option<String>,
    
    /// Optional pre-trained splice effect prediction model
    #[arg(long)]
    splice_model: Option<String>,
    
    /// Allele frequency cutoff for rare variants
    #[arg(long, default_value_t = 0.001)]
    rare_cutoff: f64,
    
    /// Output file path (supports .csv, .parquet, .json formats)
    #[arg(short, long, default_value = "annotated_variants.parquet")]
    output: String,
    
    /// Number of threads to use (0 means use all available)
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
    
    /// Chromosome to process (if omitted, process all)
    #[arg(long)]
    chromosome: Option<String>,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
    
    /// Window size for sequence context (must be odd number)
    #[arg(long, default_value_t = 5000)]
    context_size: usize,
    
    /// Export prediction confidence scores
    #[arg(long)]
    export_scores: bool,
}

/// Represents a gene interval for the Lapper interval tree
type GeneIv = Interval<GeneInfo>;

/// Gene information stored in the interval tree
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct GeneInfo {
    gene_name: String,
    gene_id: String,
    strand: String,
    biotype: String,
}

/// Genomic sequence cache to minimize reference lookups
struct SequenceCache {
    fasta_reader: Option<IndexedReader<File>>,
    cache: HashMap<String, Vec<u8>>,
    max_cache_size: usize,
}

impl SequenceCache {
    fn new(reference_path: Option<&str>, max_cache_size: usize) -> Result<Self> {
        let fasta_reader = if let Some(path) = reference_path {
            let path = Path::new(path);
            Some(
                IndexedReader::from_file(&path)
                    .with_context(|| format!("Failed to open reference genome: {}", path.display()))?
            )
        } else {
            None
        };
        
        Ok(Self {
            fasta_reader,
            cache: HashMap::new(),
            max_cache_size,
        })
    }
    
    fn fetch_sequence(&mut self, chrom: &str, pos: u64, context_size: usize) -> Result<Vec<u8>> {
        if self.fasta_reader.is_none() {
            return Err(anyhow!(AnnotationError::NoReferenceError));
        }
        
        // Calculate start and end positions for the window
        let half_size = context_size / 2;
        let start = if pos <= half_size as u64 { 0 } else { pos - half_size as u64 };
        let end = start + context_size as u64;
        
        // Create cache key
        let cache_key = format!("{}:{}-{}", chrom, start, end);
        
        // Check if sequence is in cache
        if let Some(seq) = self.cache.get(&cache_key) {
            return Ok(seq.clone());
        }
        
        // Fetch from FASTA if not in cache
        let reader = self.fasta_reader.as_mut().unwrap();
        reader.fetch(chrom, start, end)
            .with_context(|| {
                format!(
                    "Failed to fetch sequence for {}:{}-{}",
                    chrom,
                    start,
                    end
                )
            })?;
            
        // Get the sequence from the reader
        let mut sequence = Vec::new();
        reader.read(&mut sequence)
            .with_context(|| "Failed to read sequence from FASTA reader")?;
        
        // Add to cache if not too large
        if self.cache.len() < self.max_cache_size {
            self.cache.insert(cache_key, sequence.clone());
        }
        
        Ok(sequence)
    }
}

/// Variant annotation record with all computed fields
#[derive(Debug, Clone, Serialize, Deserialize)]
struct AnnotatedVariant {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    gene_name: Option<String>,
    gene_id: Option<String>,
    gene_strand: Option<String>,
    gene_biotype: Option<String>,
    gnomad_af: f64,
    is_rare: bool,
    delta_psi: Option<f64>,
    pathogenicity_score: f64,
    confidence: f64,
}

/// Build a gene interval tree from a GFF file
fn build_gene_tree<P: AsRef<Path>>(p: P) -> Result<HashMap<String, Lapper<GeneInfo>>> {
    let start_time = Instant::now();
    info!("Building gene interval trees from GFF: {:?}", p.as_ref());
    
    // Open GFF reader
    let file = File::open(&p)
        .with_context(|| format!("Failed to open GFF file: {:?}", p.as_ref()))?;
    let mut rdr = gff::reader::Reader::new(BufReader::new(file));
    
    // Create interval map per chromosome
    let mut intervals_by_chrom: HashMap<String, Vec<GeneIv>> = HashMap::new();
    
    // Process records
    let mut record_count = 0;
    let mut gene_count = 0;
    
    for record_result in rdr.records() {
        record_count += 1;
        
        // Safely unwrap record
        let record = match record_result {
            Ok(r) => r,
            Err(e) => {
                warn!("Skipping malformed GFF record: {}", e);
                continue;
            }
        };
        
        // Only process gene features
        if record.ty() != "gene" {
            continue;
        }
        
        gene_count += 1;
        
        // Extract gene information
        let gene_name = record.attributes().get("gene_name")
            .or_else(|| record.attributes().get("Name"))
            .map(|v| v.to_string())
            .unwrap_or_else(|| ".".to_string());
            
        let gene_id = record.attributes().get("gene_id")
            .or_else(|| record.attributes().get("ID"))
            .map(|v| v.to_string())
            .unwrap_or_else(|| ".".to_string());
            
        let strand = record.strand().to_string();
        let biotype = record.attributes().get("biotype")
            .or_else(|| record.attributes().get("gene_biotype"))
            .map(|v| v.to_string())
            .unwrap_or_else(|| ".".to_string());
        
        // Create gene info
        let gene_info = GeneInfo {
            gene_name,
            gene_id,
            strand,
            biotype,
        };
        
        // Create interval
        let interval = GeneIv {
            start: record.start().into(),
            stop: record.end().into(),
            val: gene_info,
        };
        
        // Add to chromosome-specific vector
        let chrom = record.reference_sequence_name().to_string();
        intervals_by_chrom.entry(chrom).or_default().push(interval);
    }
    
    // Create a Lapper for each chromosome
    let mut result = HashMap::new();
    for (chrom, intervals) in intervals_by_chrom {
        result.insert(chrom, Lapper::new(intervals));
    }
    
    let elapsed = start_time.elapsed();
    info!(
        "Built gene trees for {} chromosomes with {} genes (from {} records) in {:.2?}",
        result.len(),
        gene_count,
        record_count,
        elapsed
    );
    
    Ok(result)
}

/// Load allele frequencies from a compressed gnomAD-like file
fn load_freqs<P: AsRef<Path>>(
    bgz_path: P,
    chromosome_filter: Option<&str>,
) -> Result<FxHashMap<(String, u64, String), f64>> {
    let start_time = Instant::now();
    info!("Loading allele frequencies from {:?}", bgz_path.as_ref());
    
    let mut map = FxHashMap::default();
    let path = bgz_path.as_ref();
    
    // Open BGZF reader
    let rdr = bgzf::Reader::new(
        File::open(path).with_context(|| format!("Failed to open frequency file: {:?}", path))?,
    );
    
    // Create buffered reader
    let buf_reader = BufReader::new(rdr);
    
    // Setup progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    
    let mut line_count = 0;
    let mut loaded_count = 0;
    
    // Process lines
    for (i, line_result) in buf_reader.lines().enumerate() {
        // Update progress every 100k lines
        if i % 100_000 == 0 {
            pb.set_message(format!("Processed {} lines, loaded {} variants", i, loaded_count));
            pb.tick();
        }
        
        line_count += 1;
        
        // Safely unwrap line
        let line = match line_result {
            Ok(l) => l,
            Err(e) => {
                warn!("Error reading line from frequency file: {}", e);
                continue;
            }
        };
        
        // Skip header lines
        if line.starts_with('#') {
            continue;
        }
        
        // Parse fields
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 5 {
            warn!("Skipping malformed frequency data line: insufficient fields");
            continue;
        }
        
        // Extract chromosome
        let chrom = fields[0].to_string();
        
        // Apply chromosome filter if specified
        if let Some(target_chrom) = chromosome_filter {
            if chrom != target_chrom {
                continue;
            }
        }
        
        // Parse position
        let pos = match fields[1].parse::<u64>() {
            Ok(p) => p,
            Err(_) => {
                warn!("Skipping line with invalid position: {}", fields[1]);
                continue;
            }
        };
        
        // Extract allele
        let allele = fields[3].to_string();
        
        // Parse frequency
        let freq = match fields[4].parse::<f64>() {
            Ok(f) => f,
            Err(_) => {
                warn!("Skipping line with invalid frequency: {}", fields[4]);
                continue;
            }
        };
        
        // Insert into map
        map.insert((chrom, pos, allele), freq);
        loaded_count += 1;
    }
    
    pb.finish_with_message(format!(
        "Loaded {} frequency entries from {} lines",
        loaded_count, line_count
    ));
    
    let elapsed = start_time.elapsed();
    info!(
        "Loaded {} frequency entries in {:.2?}",
        map.len(),
        elapsed
    );
    
    Ok(map)
}

/// Perform one-hot encoding of DNA sequence for neural network input
fn one_hot_encode(_sequence: &[u8], _context_size: usize) -> Result<Tensor> {
    // Stub implementation that just returns a dummy Tensor
    Err(anyhow!("PyTorch functionality disabled"))
}

/// Predict splice effect using a pre-trained PyTorch model
fn predict_splice_effect(
    _model: &CModule,
    _sequence: &[u8],
    _context_size: usize,
) -> Result<f64> {
    // Stub implementation that returns an error
    Err(anyhow!("PyTorch functionality disabled"))
}

/// Save annotations to a file in the appropriate format
fn save_annotations(annotations: Vec<AnnotatedVariant>, output_path: &str) -> Result<()> {
    let path = Path::new(output_path);
    let extension = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    
    // Convert to DataFrame for easier output handling
    let mut df = DataFrame::new(vec![
        Series::new("chrom", annotations.iter().map(|a| a.chrom.clone()).collect::<Vec<_>>()),
        Series::new("pos", annotations.iter().map(|a| a.pos).collect::<Vec<_>>()),
        Series::new("ref_allele", annotations.iter().map(|a| a.ref_allele.clone()).collect::<Vec<_>>()),
        Series::new("alt_allele", annotations.iter().map(|a| a.alt_allele.clone()).collect::<Vec<_>>()),
        Series::new(
            "gene_name",
            annotations.iter()
                .map(|a| a.gene_name.clone().unwrap_or_else(|| "NA".to_string()))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "gene_id",
            annotations.iter()
                .map(|a| a.gene_id.clone().unwrap_or_else(|| "NA".to_string()))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "gene_strand",
            annotations.iter()
                .map(|a| a.gene_strand.clone().unwrap_or_else(|| ".".to_string()))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "gene_biotype",
            annotations.iter()
                .map(|a| a.gene_biotype.clone().unwrap_or_else(|| "NA".to_string()))
                .collect::<Vec<_>>(),
        ),
        Series::new("gnomAD_AF", annotations.iter().map(|a| a.gnomad_af).collect::<Vec<_>>()),
        Series::new("is_rare", annotations.iter().map(|a| a.is_rare).collect::<Vec<_>>()),
        Series::new(
            "delta_psi",
            annotations.iter()
                .map(|a| a.delta_psi.unwrap_or(f64::NAN))
                .collect::<Vec<_>>(),
        ),
        Series::new("pathogenicity", annotations.iter().map(|a| a.pathogenicity_score).collect::<Vec<_>>()),
        Series::new("confidence", annotations.iter().map(|a| a.confidence).collect::<Vec<_>>()),
    ])?;
    
    match extension.to_lowercase().as_str() {
        "csv" => {
            let mut file = BufWriter::new(File::create(path)?);
            CsvWriter::new(&mut file)
                .has_header(true)
                .with_delimiter(b',')
                .finish(&mut df)?;
        }
        "parquet" => {
            let mut file = File::create(path)?;
            ParquetWriter::new(&mut file)
                .with_compression(ParquetCompression::Snappy)
                .finish(&mut df)?;
        }
        "json" => {
            let json = serde_json::to_string_pretty(&annotations)?;
            let mut file = BufWriter::new(File::create(path)?);
            file.write_all(json.as_bytes())?;
        }
        _ => {
            // Default to parquet if extension not recognized
            warn!("Unrecognized file extension: {}, defaulting to Parquet format", extension);
            let output_path = format!("{}.parquet", output_path);
            let mut file = File::create(output_path)?;
            ParquetWriter::new(&mut file)
                .with_compression(ParquetCompression::Snappy)
                .finish(&mut df)?;
        }
    }
    
    info!("Saved {} annotations to {}", annotations.len(), output_path);
    Ok(())
}

fn main() -> Result<()> {
    // Measure execution time
    let start_time = Instant::now();
    
    // Parse command line arguments
    let args = Args::parse();
    
    // Configure logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(
        if args.verbose { "debug" } else { "info" },
    ))
    .format_timestamp_millis()
    .init();
    
    info!("Starting variant annotation pipeline");
    
    // Configure thread pool for parallel processing
    let num_threads = if args.threads == 0 {
        num_cpus::get()
    } else {
        args.threads
    };
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    info!("Using {} threads for parallel processing", num_threads);
    
    // Validate context size
    if args.context_size % 2 == 0 {
        return Err(anyhow!("Context size must be an odd number"));
    }
    
    // Build gene interval tree from GFF
    let gene_trees = build_gene_tree(&args.gff)?;
    
    // Load allele frequencies from gnomAD
    let freqs = load_freqs(&args.gnomad, args.chromosome.as_deref())?;
    
    // Initialize sequence cache if we have a reference
    let seq_cache = if args.splice_model.is_some() && args.reference.is_none() {
        return Err(anyhow!(AnnotationError::NoReferenceError));
    } else {
        Arc::new(Mutex::new(SequenceCache::new(args.reference.as_deref(), 1000)?))
    };
    
    // Load splice prediction model if provided - temporarily disabled
    let splice_net = None;
    if args.splice_model.is_some() {
        info!("Note: Splice model was specified, but PyTorch functionality is disabled in this build");
    }
    
    // Open VCF file
    info!("Processing variants from {}", args.vcf);
    let file = File::open(&args.vcf)
        .with_context(|| format!("Failed to open VCF file: {}", args.vcf))?;
    let mut vcf_rdr = vcf::reader::Reader::new(BufReader::new(file));
    
    // Read VCF header
    let header = vcf_rdr
        .read_header()
        .context("Failed to read VCF header")?;
    
    // Create progress bar
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    progress_bar.enable_steady_tick(std::time::Duration::from_millis(100));
    
    // Track statistics
    let stats = Arc::new(Mutex::new(HashMap::new()));
    let processed_counter = Arc::new(Mutex::new(0usize));
    
    // Process VCF records in parallel
    info!("Starting variant annotation");
    let annotations: Vec<_> = vcf_rdr
        .records(&header)
        .par_bridge()
        .filter_map(|record_result| {
            // Safely unwrap record
            let record = match record_result {
                Ok(rec) => rec,
                Err(e) => {
                    error!("Error reading VCF record: {}", e);
                    return None;
                }
            };
            
            // Extract basic variant information
            let chrom = record.chromosome().to_string();
            
            // Apply chromosome filter if specified
            if let Some(ref target_chrom) = args.chromosome {
                if chrom != *target_chrom {
                    return None;
                }
            }
            
            let pos: usize = record.position().into();
            
            // Get alleles
            let ref_allele = record.reference_bases().to_string();
            let alt_allele = record.alternate_bases().first()
                .map(|a| a.to_string())
                .unwrap_or_else(|| ".".to_string());
            
            // Lookup gene information
            let gene_info = gene_trees
                .get(&chrom)
                .and_then(|tree| tree.find(pos, pos).next())
                .map(|iv| iv.val.clone());
            
            // Get allele frequency
            let af = freqs
                .get(&(chrom.clone(), pos as u64, alt_allele.clone()))
                .copied()
                .unwrap_or(0.0);
            
            let is_rare = af < args.rare_cutoff;
            
            // Predict splice effect if model is available
            let dpsi = if splice_net.is_some() {
                // Thread-safe access to sequence cache
                let sequence = match seq_cache.lock().unwrap().fetch_sequence(&chrom, pos as u64, args.context_size) {
                    Ok(seq) => seq,
                    Err(e) => {
                        warn!("Error fetching sequence for {}:{}: {}", chrom, pos, e);
                        return None;
                    }
                };
                
                // Predict splice effect
                match predict_splice_effect(splice_net.as_ref().unwrap(), &sequence, args.context_size) {
                    Ok(effect) => Some(effect),
                    Err(e) => {
                        warn!("Error predicting splice effect for {}:{}: {}", chrom, pos, e);
                        None
                    }
                }
            } else {
                None
            };
            
            // Calculate pathogenicity score using logistic function
            // Factors: splice effect and rarity
            let dpsi_factor = dpsi.unwrap_or(0.0) * 4.0; // Scale splice effect
            let rare_factor = if is_rare { 1.0 } else { 0.0 };
            
            // Combined score through sigmoid function
            let path_score = 1.0 / (1.0 + (-dpsi_factor - rare_factor).exp());
            
            // Calculate confidence based on available data
            let confidence = if dpsi.is_some() && af > 0.0 {
                0.9 // High confidence when we have both splice prediction and frequency data
            } else if dpsi.is_some() || af > 0.0 {
                0.7 // Medium confidence with either splice prediction or frequency data
            } else {
                0.5 // Low confidence with neither
            };
            
            // Update progress and statistics
            {
                let mut count = processed_counter.lock().unwrap();
                *count += 1;
                
                if *count % 1000 == 0 {
                    progress_bar.set_message(format!("Processed {} variants", *count));
                }
                
                // Update statistics
                let mut stats_guard = stats.lock().unwrap();
                let counter = stats_guard.entry(chrom.clone()).or_insert(0);
                *counter += 1;
            }
            
            // Create annotation record
            Some(AnnotatedVariant {
                chrom,
                pos: pos as u64,
                ref_allele,
                alt_allele,
                gene_name: gene_info.as_ref().map(|g| g.gene_name.clone()),
                gene_id: gene_info.as_ref().map(|g| g.gene_id.clone()),
                gene_strand: gene_info.as_ref().map(|g| g.strand.clone()),
                gene_biotype: gene_info.as_ref().map(|g| g.biotype.clone()),
                gnomad_af: af,
                is_rare,
                delta_psi: dpsi,
                pathogenicity_score: path_score,
                confidence,
            })
        })
        .collect();
    
    // Finish progress
    progress_bar.finish_with_message(format!("Annotated {} variants", annotations.len()));
    
    // Print chromosome statistics
    info!("Annotation statistics by chromosome:");
    let stats_guard = stats.lock().unwrap();
    for (chrom, count) in stats_guard.iter() {
        info!("  {}: {} variants", chrom, count);
    }
    
    // Save annotations
    save_annotations(annotations.clone(), &args.output)?;
    
    // Print results preview
    let preview_count = std::cmp::min(annotations.len(), 12);
    if preview_count > 0 {
        println!("\nAnnotation Results Preview:");
        println!("{:-<80}", "");
        println!(
            "{:10} {:10} {:6} {:6} {:15} {:10} {:10} {:10}",
            "CHROM", "POS", "REF", "ALT", "GENE", "AF", "ΔPSI", "PATH SCORE"
        );
        println!("{:-<80}", "");
        
        for ann in annotations.iter().take(preview_count) {
            println!(
                "{:10} {:10} {:6} {:6} {:15} {:.6} {:10.3} {:10.3}",
                ann.chrom,
                ann.pos,
                ann.ref_allele,
                ann.alt_allele,
                ann.gene_name.clone().unwrap_or_default(),
                ann.gnomad_af,
                ann.delta_psi.unwrap_or(0.0),
                ann.pathogenicity_score
            );
        }
        println!("{:-<80}", "");
    }
    
    // Print timing information
    let elapsed = start_time.elapsed();
    info!(
        "Completed variant annotation in {:.2?}",
        elapsed
    );
    
    Ok(())
}