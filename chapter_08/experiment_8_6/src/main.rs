use anyhow::{Context, Result, anyhow};
use clap::{Parser, Subcommand, ValueEnum};
use ndarray::Array2;
use odgi::graph::Graph;
use onnxruntime::{
    Environment, 
    ndarray_tensor::NdArrayTensor,
    Session
};
use polars::prelude::*;
use rayon::prelude::*;
use rust_htslib::{bcf, bcf::Read};
use serde::{Serialize, Deserialize};
use whatshap_rs::phase_block;
use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::Path,
    sync::{Arc, Mutex, atomic::{AtomicUsize, Ordering}},
    time::Instant,
    thread,
};
use tracing::{info, warn, error, debug, Level};
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use thiserror::Error;
use tempfile::NamedTempFile;

/// Custom error types for the variant scoring pipeline
#[derive(Error, Debug)]
pub enum ScoringError {
    #[error("Failed to load graph: {0}")]
    GraphLoadError(String),
    
    #[error("Failed to load model: {0}")]
    ModelLoadError(String),
    
    #[error("Failed to process VCF: {0}")]
    VcfProcessError(String),
    
    #[error("Failed to run inference: {0}")]
    InferenceError(String),
    
    #[error("Failed to phase variants: {0}")]
    PhasingError(String),
    
    #[error("Failed to write output: {0}")]
    OutputError(String),
    
    #[error("Invalid input data: {0}")]
    InvalidInputError(String),
}

/// Supported output formats
#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Parquet,
    Ipc,
    Csv,
    Json,
    Tsv,
}

/// Command line interface
#[derive(Parser, Debug)]
#[clap(
    name = "variant-scorer",
    about = "Score genetic variants using pangenome graphs and machine learning",
    version
)]
struct Cli {
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
    
    /// Number of threads to use (0 = use all available)
    #[arg(short, long, default_value = "0")]
    threads: usize,
    
    #[command(subcommand)]
    cmd: Command,
}

/// Available commands
#[derive(Subcommand, Debug)]
enum Command {
    /// Score variants using a pangenome graph and ML model
    Score {
        /// Path to pangenome graph in ODGI format
        #[arg(long)]
        graph: String,
        
        /// Path to VCF file with variants to score
        #[arg(long)]
        vcf: String,
        
        /// Path to ONNX model for scoring
        #[arg(long)]
        model: String,
        
        /// Path to output file
        #[arg(long)]
        out: String,
        
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Ipc)]
        format: OutputFormat,
        
        /// Batch size for processing
        #[arg(long, default_value = "1000")]
        batch_size: usize,
        
        /// Window size for phasing (in bp)
        #[arg(long, default_value = "1000")]
        phase_window: i32,
        
        /// Skip phasing step
        #[arg(long)]
        skip_phasing: bool,
        
        /// Include additional features from graph
        #[arg(long)]
        extended_features: bool,
        
        /// Filter out variants with score below threshold
        #[arg(long)]
        min_score: Option<f32>,
    },
    
    /// Batch score variants from multiple VCFs
    BatchScore {
        /// Path to pangenome graph in ODGI format
        #[arg(long)]
        graph: String,
        
        /// Path to file containing list of VCF files
        #[arg(long)]
        vcf_list: String,
        
        /// Path to ONNX model for scoring
        #[arg(long)]
        model: String,
        
        /// Output directory
        #[arg(long)]
        out_dir: String,
        
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Ipc)]
        format: OutputFormat,
        
        /// Window size for phasing (in bp)
        #[arg(long, default_value = "1000")]
        phase_window: i32,
    },
}

/// Configuration for scoring process
#[derive(Debug, Clone)]
struct ScoringConfig {
    batch_size: usize,
    phase_window: i32,
    skip_phasing: bool,
    extended_features: bool,
    min_score: Option<f32>,
    output_format: OutputFormat,
}

/// Default configuration
impl Default for ScoringConfig {
    fn default() -> Self {
        Self {
            batch_size: 1000,
            phase_window: 1000,
            skip_phasing: false,
            extended_features: false,
            min_score: None,
            output_format: OutputFormat::Ipc,
        }
    }
}

/// Variant information for scoring
#[derive(Debug, Clone, Serialize, Deserialize)]
struct VariantInfo {
    chrom: String,
    pos: i64,
    ref_allele: String,
    alt_allele: String,
    score: f64,
    phase_block: String,
    node_id: Option<u64>,
    node_degree: Option<u32>,
    centrality: Option<f64>,
}

/// Statistics for reporting
#[derive(Debug, Default, Serialize, Deserialize)]
struct ScoringStats {
    total_variants: usize,
    processed_variants: usize,
    filtered_variants: usize,
    high_scoring_variants: usize,
    multi_allelic_variants: usize,
    phased_variants: usize,
    elapsed_seconds: f64,
}

/// Main entry point
fn main() -> Result<()> {
    // Parse command line arguments
    let cli = Cli::parse();
    
    // Configure logging
    let log_level = if cli.verbose { Level::DEBUG } else { Level::INFO };
    tracing_subscriber::fmt()
        .with_max_level(log_level)
        .init();
    
    // Configure thread pool
    let num_threads = if cli.threads == 0 {
        num_cpus::get()
    } else {
        cli.threads
    };
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    info!("Using {} threads for parallel processing", num_threads);
    
    // Start timing
    let start_time = Instant::now();
    
    // Execute command
    let result = match &cli.cmd {
        Command::Score {
            graph,
            vcf,
            model,
            out,
            format,
            batch_size,
            phase_window,
            skip_phasing,
            extended_features,
            min_score,
        } => {
            let config = ScoringConfig {
                batch_size: *batch_size,
                phase_window: *phase_window,
                skip_phasing: *skip_phasing,
                extended_features: *extended_features,
                min_score: *min_score,
                output_format: *format,
            };
            
            run_score(graph, vcf, model, out, &config)
        }
        
        Command::BatchScore {
            graph,
            vcf_list,
            model,
            out_dir,
            format,
            phase_window,
        } => {
            let config = ScoringConfig {
                batch_size: 1000,
                phase_window: *phase_window,
                skip_phasing: false,
                extended_features: true,
                min_score: None,
                output_format: *format,
            };
            
            run_batch_score(graph, vcf_list, model, out_dir, &config)
        }
    };
    
    // Log execution time
    let elapsed = start_time.elapsed();
    info!("Total execution time: {:.2?}", elapsed);
    
    result
}

/// Load and validate an ODGI pangenome graph
fn load_graph(graph_path: &str) -> Result<Graph> {
    info!("Loading pangenome graph from: {}", graph_path);
    let start = Instant::now();
    
    let graph = Graph::from_json_path(graph_path)
        .map_err(|e| anyhow!(ScoringError::GraphLoadError(e)))?;
    
    let node_count = graph.node_count();
    let edge_count = graph.edge_count();
    
    info!(
        "Loaded graph with {} nodes and {} edges in {:.2?}",
        node_count,
        edge_count,
        start.elapsed()
    );
    
    // Validate graph has content
    if node_count == 0 {
        return Err(anyhow!(ScoringError::GraphLoadError(
            "Graph contains no nodes".to_string()
        )));
    }
    
    Ok(graph)
}

/// Initialize ONNX runtime and load model
fn load_model(model_path: &str) -> Result<(Environment, Session)> {
    info!("Loading ONNX model from: {}", model_path);
    let start = Instant::now();
    
    // Initialize ONNX runtime environment
    let environment = Environment::builder()
        .with_name("variant_scorer")
        .build()
        .context("Failed to build ONNX environment")?;
    
    // Create session with optimized execution providers
    let mut session_builder = environment.new_session_builder()?;
    
    // Check for GPU availability and configure execution providers
    #[cfg(feature = "cuda")]
    {
        let cuda_provider = onnxruntime::ExecutionProvider::CUDA(Default::default());
        session_builder = session_builder.with_execution_providers([cuda_provider])?;
        debug!("Using CUDA execution provider for ONNX inference");
    }
    
    #[cfg(not(feature = "cuda"))]
    {
        let cpu_provider = onnxruntime::ExecutionProvider::CPU(Default::default());
        session_builder = session_builder.with_execution_providers([cpu_provider])?;
        debug!("Using CPU execution provider for ONNX inference");
    }
    
    // Load the model
    let session = session_builder
        .with_model_from_file(model_path)
        .with_context(|| format!("Failed to load ONNX model from {}", model_path))?;
    
    // Get model metadata
    let model_metadata = session.model_metadata()?;
    let input_names = model_metadata.inputs.iter().map(|i| i.name.clone()).collect::<Vec<_>>();
    let output_names = model_metadata.outputs.iter().map(|o| o.name.clone()).collect::<Vec<_>>();
    
    info!(
        "Loaded ONNX model in {:.2?} with inputs: {:?}, outputs: {:?}",
        start.elapsed(),
        input_names,
        output_names
    );
    
    Ok((environment, session))
}

/// Run inference on a batch of variants
fn run_inference(
    session: &Session,
    features: Array2<f32>,
    extended_features: bool,
) -> Result<Vec<f32>> {
    // Validate feature array dimensions
    let expected_features = if extended_features { 5 } else { 3 };
    if features.shape()[1] != expected_features {
        return Err(anyhow!(ScoringError::InferenceError(format!(
            "Invalid feature dimensions: expected {} features, got {}",
            expected_features,
            features.shape()[1]
        ))));
    }
    
    // Create input tensor
    let input_tensor = NdArrayTensor::from_array(features);
    
    // Run inference
    let outputs = session
        .run(vec![input_tensor])
        .context("Failed to run ONNX inference")?;
    
    // Extract scores from output tensor
    let scores: Vec<f32> = outputs[0]
        .float_array()
        .context("Failed to get float array from ONNX output")?
        .iter()
        .copied()
        .collect();
    
    Ok(scores)
}

/// Extract features from a variant and graph context
fn extract_features(
    graph: &Graph,
    chrom: &str,
    pos: i64,
    ref_allele: &str,
    alt_allele: &str,
    extended_features: bool,
) -> Result<Vec<f32>> {
    // Basic features: reference length, alternate length
    let mut features = vec![
        ref_allele.len() as f32,
        alt_allele.len() as f32,
    ];
    
    // Get graph context at this position
    let node_degree = graph.degree_at(chrom, pos as u64).unwrap_or(0) as f32;
    features.push(node_degree);
    
    // Add extended features if requested
    if extended_features {
        // Get node centrality (proxy for importance in graph)
        let centrality = graph.centrality_at(chrom, pos as u64).unwrap_or(0.0) as f32;
        features.push(centrality);
        
        // Compute sequence complexity feature
        // Simple implementation: ratio of unique k-mers to length
        let seq_complexity = compute_sequence_complexity(alt_allele);
        features.push(seq_complexity);
    }
    
    Ok(features)
}

/// Compute sequence complexity (simple k-mer based approach)
fn compute_sequence_complexity(sequence: &str) -> f32 {
    if sequence.len() <= 3 {
        return 1.0;
    }
    
    let k = 3; // k-mer size
    let mut kmers = std::collections::HashSet::new();
    
    for i in 0..=(sequence.len() - k) {
        kmers.insert(&sequence[i..(i + k)]);
    }
    
    // Ratio of unique k-mers to possible k-mers
    let max_kmers = sequence.len() - k + 1;
    kmers.len() as f32 / max_kmers as f32
}

/// Score variants in a VCF file
fn run_score(
    graph_path: &str,
    vcf_path: &str,
    model_path: &str,
    out_path: &str,
    config: &ScoringConfig,
) -> Result<()> {
    let start_time = Instant::now();
    
    // Load graph
    let graph = load_graph(graph_path)?;
    
    // Load model
    let (_environment, session) = load_model(model_path)?;
    
    // Setup progress tracking
    let multi_progress = MultiProgress::new();
    let main_progress = multi_progress.add(ProgressBar::new_spinner());
    main_progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    
    // Spawn a thread to render progress bars
    let progress_thread = thread::spawn(move || {
        // In our mock implementation, we don't need to call join()
        // Just sleep for a bit to let the progress bars update
        thread::sleep(std::time::Duration::from_millis(100));
    });
    
    main_progress.set_message(format!("Processing variants from {}", vcf_path));
    
    // Validate VCF exists
    if !Path::new(vcf_path).exists() {
        return Err(anyhow!(ScoringError::VcfProcessError(format!(
            "VCF file not found: {}",
            vcf_path
        ))));
    }
    
    // Open VCF reader
    let mut reader = bcf::Reader::from_path(vcf_path)
        .with_context(|| format!("Failed to open VCF file: {}", vcf_path))?;
    
    // Count total variants for progress tracking
    let total_variants = count_variants(vcf_path)?;
    let batch_progress = multi_progress.add(
        ProgressBar::new(total_variants as u64)
            .with_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} variants ({eta})")
                    .unwrap()
                    .progress_chars("=>-"),
            ),
    );
    
    // Setup statistics tracking
    let stats = Arc::new(Mutex::new(ScoringStats {
        total_variants,
        ..Default::default()
    }));
    
    // Setup shared data structure for collecting results
    let variants_info = Arc::new(Mutex::new(Vec::with_capacity(total_variants)));
    let counter = Arc::new(AtomicUsize::new(0));
    
    // Process VCF in batches
    let batch_size = config.batch_size;
    let mut batch = Vec::with_capacity(batch_size);
    let mut records = reader.records();
    
    // Get batch of records
    while let Some(record_result) = records.next() {
        let record = record_result.with_context(|| "Failed to read VCF record")?;
        batch.push(record);
        
        if batch.len() >= batch_size {
            // Process batch
            process_batch(
                &batch,
                &graph,
                &session,
                config,
                &variants_info,
                &stats,
                &counter,
                &batch_progress,
            )?;
            
            // Clear batch
            batch.clear();
        }
    }
    
    // Process final batch if there are remaining records
    if !batch.is_empty() {
        process_batch(
            &batch,
            &graph,
            &session,
            config,
            &variants_info,
            &stats,
            &counter,
            &batch_progress,
        )?;
    }
    
    // Get final count
    let processed_count = counter.load(Ordering::SeqCst);
    batch_progress.finish_with_message(format!("Processed {} variants", processed_count));
    
    // Update elapsed time in stats
    {
        let mut stats_guard = stats.lock().unwrap();
        stats_guard.elapsed_seconds = start_time.elapsed().as_secs_f64();
        stats_guard.processed_variants = processed_count;
    }
    
    // Get all variant information
    let all_variants = {
        let guard = variants_info.lock().unwrap();
        guard.clone()
    };
    
    // Save results
    main_progress.set_message(format!("Writing results to {}", out_path));
    save_results(&all_variants, out_path, config.output_format)?;
    
    // Print statistics
    let stats_guard = stats.lock().unwrap();
    print_statistics(&stats_guard);
    
    // Finish progress
    main_progress.finish_with_message(format!(
        "Completed scoring {} variants in {:.2?}",
        processed_count,
        start_time.elapsed()
    ));
    
    // Wait for progress thread to finish
    drop(batch_progress);
    progress_thread.join().unwrap();
    
    Ok(())
}

/// Count variants in a VCF file
fn count_variants(vcf_path: &str) -> Result<usize> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let count = reader.records().count();
    Ok(count)
}

/// Process a batch of variants
fn process_batch(
    batch: &[bcf::Record],
    graph: &Graph,
    session: &Session,
    config: &ScoringConfig,
    variants_info: &Arc<Mutex<Vec<VariantInfo>>>,
    stats: &Arc<Mutex<ScoringStats>>,
    counter: &Arc<AtomicUsize>,
    progress: &ProgressBar,
) -> Result<()> {
    // Create batch feature matrix
    let mut feature_vectors = Vec::with_capacity(batch.len());
    let mut variant_meta = Vec::with_capacity(batch.len());
    
    // Process each variant in the batch
    for record in batch {
        // Get chromosome and position
        // Use rid to get chromosome name since chrom() method doesn't exist
        let rid = record.rid().ok_or_else(|| anyhow!("Record has no RID"))?;
        let header = record.header();
        let chrom = std::str::from_utf8(header.rid2name(rid)?)
            .context("Failed to decode chromosome name")?
            .to_owned();
        let pos = record.pos();
        
        // Get alleles - accessing directly without error matching since it returns Vec<&[u8]>
        let alleles = record.alleles();
        
        // Skip if not biallelic
        if alleles.len() != 2 {
            // Update multi-allelic counter
            {
                let mut stats_guard = stats.lock().unwrap();
                stats_guard.multi_allelic_variants += 1;
            }
            continue;
        }
        
        // Convert alleles to strings
        let ref_allele = std::str::from_utf8(alleles[0])
            .context("Failed to decode reference allele")?
            .to_owned();
        let alt_allele = std::str::from_utf8(alleles[1])
            .context("Failed to decode alternate allele")?
            .to_owned();
        
        // Extract features
        match extract_features(
            graph,
            &chrom,
            pos as i64,
            &ref_allele,
            &alt_allele,
            config.extended_features,
        ) {
            Ok(features) => {
                // Store features and metadata
                feature_vectors.push(features);
                variant_meta.push((chrom, pos as i64, ref_allele, alt_allele));
            }
            Err(err) => {
                warn!("Failed to extract features for variant at {}:{}: {}", chrom, pos, err);
                continue;
            }
        }
    }
    
    // Skip if no valid variants
    if feature_vectors.is_empty() {
        return Ok(());
    }
    
    // Create feature array
    let feature_dim = if config.extended_features { 5 } else { 3 };
    let mut feature_array = Array2::zeros((feature_vectors.len(), feature_dim));
    
    for (i, features) in feature_vectors.iter().enumerate() {
        for (j, &value) in features.iter().enumerate() {
            feature_array[[i, j]] = value;
        }
    }
    
    // Run inference
    let scores = run_inference(session, feature_array, config.extended_features)?;
    
    // Phase variants if requested
    let phase_results = if !config.skip_phasing {
        // In our simplified implementation, we'll just phase each variant directly
        batch
            .par_iter()
            .map(|record| {
                let mut record_copy = record.clone();
                match phase_block(&mut record_copy, config.phase_window) {
                    Ok(phase_tag) => {
                        // Increment phased counter if tag is not empty
                        if !phase_tag.is_empty() && phase_tag != "." {
                            let mut stats_guard = stats.lock().unwrap();
                            stats_guard.phased_variants += 1;
                        }
                        Ok(phase_tag)
                    }
                    Err(err) => Err(anyhow!(ScoringError::PhasingError(format!(
                        "Failed to phase variant: {}",
                        err
                    )))),
                }
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        // If phasing is skipped, just use empty tags
        vec![".".to_string(); batch.len()]
    };
    
    // Create variant info records
    let mut new_variants = Vec::with_capacity(feature_vectors.len());
    
    for (i, (chrom, pos, ref_allele, alt_allele)) in variant_meta.into_iter().enumerate() {
        let score = scores[i];
        
        // Skip if below threshold
        if let Some(min_score) = config.min_score {
            if score < min_score {
                // Update filtered counter
                {
                    let mut stats_guard = stats.lock().unwrap();
                    stats_guard.filtered_variants += 1;
                }
                continue;
            }
        }
        
        // Get node ID and additional graph features
        let node_id = graph.node_at(&chrom, pos as u64);
        let node_degree = node_id.map(|id| graph.degree(id));
        let centrality = node_id.map(|id| graph.centrality(id));
        
        // Update high scoring counter
        if score >= 0.7 {
            let mut stats_guard = stats.lock().unwrap();
            stats_guard.high_scoring_variants += 1;
        }
        
        // Create variant info
        let variant_info = VariantInfo {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            score: score as f64,
            phase_block: phase_results[i].clone(),
            node_id,
            node_degree,
            centrality,
        };
        
        new_variants.push(variant_info);
    }
    
    // Update progress
    let new_count = new_variants.len();
    counter.fetch_add(new_count, Ordering::SeqCst);
    progress.inc(new_count as u64);
    
    // Add variants to the shared collection
    {
        let mut variants_guard = variants_info.lock().unwrap();
        variants_guard.extend(new_variants);
    }
    
    Ok(())
}

/// Run batch scoring on multiple VCF files
fn run_batch_score(
    graph_path: &str,
    vcf_list_path: &str,
    model_path: &str,
    out_dir: &str,
    config: &ScoringConfig,
) -> Result<()> {
    let start_time = Instant::now();
    
    // Load graph
    let _graph = load_graph(graph_path)?;
    
    // Load model
    let (_environment, _session) = load_model(model_path)?;
    
    // Read VCF list
    let vcf_files = read_file_list(vcf_list_path)
        .with_context(|| format!("Failed to read VCF list from {}", vcf_list_path))?;
    
    // Create output directory if it doesn't exist
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("Failed to create output directory: {}", out_dir))?;
    
    // Process each VCF file
    for (idx, vcf_path) in vcf_files.iter().enumerate() {
        info!("Processing file {}/{}: {}", idx + 1, vcf_files.len(), vcf_path);
        
        // Create output path
        let file_name = Path::new(vcf_path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .replace(".vcf", "")
            .replace(".gz", "");
        
        let out_path = format!(
            "{}/{}.scored.{}",
            out_dir,
            file_name,
            match config.output_format {
                OutputFormat::Parquet => "parquet",
                OutputFormat::Ipc => "arrow",
                OutputFormat::Csv => "csv",
                OutputFormat::Json => "json",
                OutputFormat::Tsv => "tsv",
            }
        );
        
        // Process this VCF
        match run_score(graph_path, vcf_path, model_path, &out_path, config) {
            Ok(_) => info!("Successfully processed {}", vcf_path),
            Err(e) => {
                error!("Failed to process {}: {}", vcf_path, e);
                // Continue with next file
            }
        }
    }
    
    info!(
        "Batch processing completed in {:.2?}",
        start_time.elapsed()
    );
    
    Ok(())
}

/// Read a list of files from a text file
fn read_file_list(path: &str) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("Failed to open file list: {}", path))?;
    let reader = BufReader::new(file);
    let mut file_list = Vec::new();
    
    for line in io::BufRead::lines(reader) {
        let line = line?;
        let trimmed = line.trim();
        
        // Skip empty lines and comments
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            file_list.push(trimmed.to_string());
        }
    }
    
    Ok(file_list)
}

/// Save results to a file in the specified format
fn save_results(
    variants: &[VariantInfo],
    out_path: &str,
    format: OutputFormat,
) -> Result<()> {
    // Create a temporary file for writing
    let dir = Path::new(out_path).parent().unwrap_or_else(|| Path::new("."));
    let temp_file = NamedTempFile::new_in(dir)?;
    
    // Get a mutable reference to the DataFrame before writing
    let mut df = create_dataframe(variants)?;
    
    match format {
        OutputFormat::Parquet => {
            ParquetWriter::new(File::create(temp_file.path())?)
                .with_compression(ParquetCompression::Snappy)
                .finish(&mut df)?;
        }
        OutputFormat::Ipc => {
            IpcWriter::new(File::create(temp_file.path())?)
                .finish(&mut df)?;
        }
        OutputFormat::Csv => {
            CsvWriter::new(File::create(temp_file.path())?)
                .has_header(true)
                .with_delimiter(b',')
                .finish(&mut df)?;
        }
        OutputFormat::Json => {
            // Serialize to JSON using serde
            let json = serde_json::to_string_pretty(variants)?;
            let mut file = BufWriter::new(File::create(temp_file.path())?);
            file.write_all(json.as_bytes())?;
            file.flush()?;
        }
        OutputFormat::Tsv => {
            CsvWriter::new(File::create(temp_file.path())?)
                .has_header(true)
                .with_delimiter(b'\t')
                .finish(&mut df)?;
        }
    }
    
    // Rename temporary file to the target path (atomic operation)
    temp_file.persist(out_path)
        .with_context(|| format!("Failed to write output file: {}", out_path))?;
    
    info!("Results saved to {}", out_path);
    Ok(())
}

/// Create a DataFrame from variant information
fn create_dataframe(variants: &[VariantInfo]) -> Result<DataFrame> {
    // Create vectors for each column
    let chroms = variants.iter().map(|v| v.chrom.clone()).collect::<Vec<_>>();
    let positions = variants.iter().map(|v| v.pos).collect::<Vec<_>>();
    let ref_alleles = variants.iter().map(|v| v.ref_allele.clone()).collect::<Vec<_>>();
    let alt_alleles = variants.iter().map(|v| v.alt_allele.clone()).collect::<Vec<_>>();
    let scores = variants.iter().map(|v| v.score).collect::<Vec<_>>();
    let phase_blocks = variants.iter().map(|v| v.phase_block.clone()).collect::<Vec<_>>();
    
    // Optional columns
    let node_ids = variants.iter().map(|v| v.node_id.unwrap_or(0)).collect::<Vec<_>>();
    let node_degrees = variants.iter().map(|v| v.node_degree.unwrap_or(0)).collect::<Vec<_>>();
    let centralities = variants.iter().map(|v| v.centrality.unwrap_or(0.0)).collect::<Vec<_>>();
    
    // Create DataFrame
    let df_columns = vec![
        Series::new("chrom", chroms),
        Series::new("pos", positions),
        Series::new("ref", ref_alleles),
        Series::new("alt", alt_alleles),
        Series::new("score", scores),
        Series::new("phase", phase_blocks),
        Series::new("node_id", node_ids),
        Series::new("node_degree", node_degrees),
        Series::new("centrality", centralities),
    ];
    
    DataFrame::new(df_columns)
        .with_context(|| "Failed to create DataFrame from variant data")
}

/// Print statistics about the scoring process
fn print_statistics(stats: &ScoringStats) {
    println!("\n===== Variant Scoring Statistics =====");
    println!("Total variants: {}", stats.total_variants);
    println!("Processed variants: {}", stats.processed_variants);
    println!("High scoring variants (≥0.7): {}", stats.high_scoring_variants);
    println!("Filtered variants: {}", stats.filtered_variants);
    println!("Multi-allelic variants: {}", stats.multi_allelic_variants);
    println!("Phased variants: {}", stats.phased_variants);
    println!("Processing time: {:.2} seconds", stats.elapsed_seconds);
    println!("=====================================\n");
}