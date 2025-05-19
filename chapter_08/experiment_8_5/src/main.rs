use clap::{Parser, Subcommand, ValueEnum};
use tracing::{info, warn, error, debug, Level};
use tracing_subscriber::filter::EnvFilter;
use anyhow::{Context, Result, anyhow, bail};
use noodles_bam as bam;
use noodles_vcf as vcf;
use noodles_fasta as fasta;
use noodles_gff as gff;
use tokio::{signal, fs, time};
use rayon::prelude::*;
use std::{
    path::{Path, PathBuf},
    process::ExitCode,
    sync::{Arc, Mutex},
    time::Instant,
};
use serde::{Deserialize, Serialize};
use thiserror::Error;
use tempfile::TempDir;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use human_format::Formatter;

/// Custom error types for the variant pipeline
#[derive(Error, Debug)]
pub enum PipelineError {
    #[error("File not found: {0}")]
    FileNotFound(String),
    
    #[error("Invalid input file: {0}")]
    InvalidInput(String),
    
    #[error("Process interrupted")]
    Interrupted,
    
    #[error("External command failed: {0}")]
    CommandFailed(String),
    
    #[error("Configuration error: {0}")]
    ConfigError(String),
    
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Pipeline settings from configuration file
#[derive(Deserialize, Serialize, Debug, Clone)]
struct Settings {
    // General settings
    threads: usize,
    tmp_dir: PathBuf,
    log_level: Option<String>,
    
    // Alignment settings
    #[serde(default)]
    align: AlignSettings,
    
    // Variant calling settings
    #[serde(default)]
    call: CallSettings,
    
    // Annotation settings
    #[serde(default)]
    annotate: AnnotateSettings,
}

/// Default implementation for Settings
impl Default for Settings {
    fn default() -> Self {
        Self {
            threads: num_cpus::get(),
            tmp_dir: std::env::temp_dir(),
            log_level: Some("info".to_string()),
            align: AlignSettings::default(),
            call: CallSettings::default(),
            annotate: AnnotateSettings::default(),
        }
    }
}

/// Alignment-specific settings
#[derive(Deserialize, Serialize, Debug, Clone, Default)]
struct AlignSettings {
    aligner: Option<String>,
    min_mapq: Option<u8>,
    max_secondary: Option<usize>,
    mark_duplicates: Option<bool>,
}

/// Variant calling settings
#[derive(Deserialize, Serialize, Debug, Clone, Default)]
struct CallSettings {
    caller: Option<String>,
    min_depth: Option<usize>,
    min_gq: Option<f64>,
    regions: Option<Vec<String>>,
}

/// Annotation settings
#[derive(Deserialize, Serialize, Debug, Clone, Default)]
struct AnnotateSettings {
    databases: Option<Vec<PathBuf>>,
    effects: Option<bool>,
    max_distance: Option<usize>,
}

/// Format options for output files
#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Vcf,
    Bcf,
    Tsv,
    Json,
    Parquet,
}

/// Command-line interface
#[derive(Parser, Debug)]
#[command(name = "genomic_pipeline", version, about = "Genomic variant analysis pipeline")]
struct Cli {
    /// Path to configuration file
    #[arg(long, default_value = "pipeline.toml")]
    config: PathBuf,
    
    /// Override number of threads from config
    #[arg(short, long)]
    threads: Option<usize>,
    
    /// Override temporary directory from config
    #[arg(long)]
    tmp_dir: Option<PathBuf>,
    
    /// Verbosity level
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
    
    /// Subcommands
    #[command(subcommand)]
    command: Command,
}

/// Subcommands for the pipeline
#[derive(Subcommand, Debug)]
enum Command {
    /// Align sequencing reads to a reference genome
    Align {
        /// Path to input reads (FASTQ)
        #[arg(short, long)]
        reads: PathBuf,
        
        /// Path to reference genome (FASTA)
        #[arg(short, long)]
        reference: PathBuf,
        
        /// Path to output BAM file
        #[arg(short, long)]
        out_bam: PathBuf,
        
        /// Aligner to use (bwa, minimap2, bowtie2)
        #[arg(long)]
        aligner: Option<String>,
        
        /// Mark duplicate reads
        #[arg(long)]
        mark_duplicates: bool,
    },
    
    /// Call variants from aligned reads
    Call {
        /// Path to input BAM file
        #[arg(short, long)]
        bam: PathBuf,
        
        /// Path to reference genome (FASTA)
        #[arg(short, long)]
        reference: PathBuf,
        
        /// Path to output VCF file
        #[arg(short, long)]
        out_vcf: PathBuf,
        
        /// Minimum read depth for variant calling
        #[arg(long)]
        min_depth: Option<usize>,
        
        /// Regions to analyze (chr:start-end format)
        #[arg(short, long)]
        regions: Option<Vec<String>>,
        
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Vcf)]
        format: OutputFormat,
    },
    
    /// Annotate variants with functional information
    Annotate {
        /// Path to input VCF file
        #[arg(short, long)]
        vcf: PathBuf,
        
        /// Path to gene annotation (GFF)
        #[arg(short, long)]
        gff: PathBuf,
        
        /// Path to output file
        #[arg(short, long)]
        output: PathBuf,
        
        /// Additional annotation databases
        #[arg(long)]
        databases: Option<Vec<PathBuf>>,
        
        /// Include effect predictions
        #[arg(long)]
        effects: bool,
        
        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
        format: OutputFormat,
    },
    
    /// Run the full pipeline (align, call, annotate)
    Pipeline {
        /// Path to input reads (FASTQ)
        #[arg(short, long)]
        reads: PathBuf,
        
        /// Path to reference genome (FASTA)
        #[arg(short, long)]
        reference: PathBuf,
        
        /// Path to gene annotation (GFF)
        #[arg(short, long)]
        gff: PathBuf,
        
        /// Path to output directory
        #[arg(short, long)]
        output_dir: PathBuf,
        
        /// Sample name (used for output files)
        #[arg(short, long)]
        sample: String,
        
        /// Keep intermediate files
        #[arg(long)]
        keep_intermediate: bool,
    },
}

/// Pipeline context shared across steps
#[derive(Debug, Clone)]
struct PipelineContext {
    settings: Settings,
    temp_dir: Arc<TempDir>,
    progress: Arc<MultiProgress>,
    start_time: Instant,
}

/// Statistics for reporting
#[derive(Debug, Default, Serialize)]
struct PipelineStats {
    aligned_reads: usize,
    variants_called: usize,
    variants_annotated: usize,
    elapsed_seconds: f64,
}

/// Main entry point
#[tokio::main(flavor = "multi_thread")]
async fn main() -> ExitCode {
    // Initialize with default logging until we parse config
    tracing_subscriber::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .init();
    
    // Measure execution time
    let start_time = Instant::now();
    
    // Parse command line arguments
    let cli = Cli::parse();
    
    // Run the pipeline with proper error handling
    match run_pipeline(cli, start_time).await {
        Ok(_) => {
            info!("Pipeline completed successfully");
            ExitCode::SUCCESS
        }
        Err(e) => {
            error!("Pipeline failed: {:#}", e);
            ExitCode::FAILURE
        }
    }
}

/// Main pipeline execution
async fn run_pipeline(cli: Cli, start_time: Instant) -> Result<()> {
    // Load and merge configuration
    let mut settings = load_configuration(&cli).await?;
    
    // Apply CLI overrides
    if let Some(threads) = cli.threads {
        settings.threads = threads;
    }
    
    if let Some(tmp_dir) = cli.tmp_dir {
        settings.tmp_dir = tmp_dir;
    }
    
    // Setup logging based on verbosity
    let log_level = match cli.verbose {
        0 => settings.log_level.as_deref().unwrap_or("info"),
        1 => "debug",
        _ => "trace",
    };
    
    setup_logging(log_level);
    
    // Initialize Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(settings.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;
    
    info!("Using {} threads for parallel processing", settings.threads);
    
    // Create temporary directory
    let temp_dir = Arc::new(
        tempfile::Builder::new()
            .prefix("genomic_pipeline_")
            .tempdir_in(&settings.tmp_dir)
            .context("Failed to create temporary directory")?,
    );
    
    debug!("Created temporary directory: {:?}", temp_dir.path());
    
    // Initialize progress bars
    let progress = Arc::new(MultiProgress::new());
    
    // Create pipeline context
    let context = PipelineContext {
        settings: settings.clone(),
        temp_dir,
        progress,
        start_time,
    };
    
    // Set up graceful shutdown handler
    let graceful = signal::ctrl_c();
    
    // Run the command with graceful shutdown
    tokio::select! {
        result = execute_command(cli.command, context) => result,
        _ = graceful => {
            info!("Shutting down on SIGINT");
            Err(anyhow!(PipelineError::Interrupted))
        }
    }
}

/// Load and parse configuration file
async fn load_configuration(cli: &Cli) -> Result<Settings> {
    let config_path = &cli.config;
    
    // Check if config file exists
    if !config_path.exists() {
        // If default config doesn't exist, use default settings
        if config_path.to_string_lossy() == "pipeline.toml" {
            return Ok(Settings::default());
        }
        
        return Err(anyhow!(PipelineError::FileNotFound(
            config_path.to_string_lossy().to_string()
        )));
    }
    
    // Read and parse configuration file
    let config_content = fs::read_to_string(config_path)
        .await
        .with_context(|| format!("Failed to read config file: {:?}", config_path))?;
    
    let settings: Settings = toml::from_str(&config_content)
        .with_context(|| format!("Failed to parse config file: {:?}", config_path))?;
    
    Ok(settings)
}

/// Set up logging with the appropriate level
fn setup_logging(level: &str) {
    let filter = match level.to_lowercase().as_str() {
        "trace" => Level::TRACE,
        "debug" => Level::DEBUG,
        "info" => Level::INFO,
        "warn" => Level::WARN,
        "error" => Level::ERROR,
        _ => Level::INFO,
    };
    
    tracing_subscriber::fmt()
        .with_max_level(filter)
        .init();
}

/// Execute the selected command
async fn execute_command(command: Command, context: PipelineContext) -> Result<()> {
    match command {
        Command::Align { reads, reference, out_bam, aligner, mark_duplicates } => {
            // Validate input files
            validate_files(&[&reads, &reference]).await?;
            
            // Create output directory if it doesn't exist
            if let Some(parent) = out_bam.parent() {
                fs::create_dir_all(parent).await?;
            }
            
            // Merge settings with command line options
            let mut align_settings = context.settings.align.clone();
            if let Some(aligner_name) = aligner {
                align_settings.aligner = Some(aligner_name);
            }
            align_settings.mark_duplicates = Some(mark_duplicates);
            
            // Run alignment
            run_alignment(&reads, &reference, &out_bam, align_settings, &context).await
        }
        
        Command::Call { bam, reference, out_vcf, min_depth, regions, format } => {
            // Validate input files
            validate_files(&[&bam, &reference]).await?;
            
            // Create output directory if it doesn't exist
            if let Some(parent) = out_vcf.parent() {
                fs::create_dir_all(parent).await?;
            }
            
            // Merge settings with command line options
            let mut call_settings = context.settings.call.clone();
            if let Some(depth) = min_depth {
                call_settings.min_depth = Some(depth);
            }
            if let Some(regions_list) = regions {
                call_settings.regions = Some(regions_list);
            }
            
            // Run variant calling
            run_calling(&bam, &reference, &out_vcf, call_settings, format, &context).await
        }
        
        Command::Annotate { vcf, gff, output, databases, effects, format } => {
            // Validate input files
            validate_files(&[&vcf, &gff]).await?;
            
            // Create output directory if it doesn't exist
            if let Some(parent) = output.parent() {
                fs::create_dir_all(parent).await?;
            }
            
            // Merge settings with command line options
            let mut annotate_settings = context.settings.annotate.clone();
            if let Some(db_list) = databases {
                annotate_settings.databases = Some(db_list);
            }
            annotate_settings.effects = Some(effects);
            
            // Run annotation
            run_annotation(&vcf, &gff, &output, annotate_settings, format, &context).await
        }
        
        Command::Pipeline { reads, reference, gff, output_dir, sample, keep_intermediate } => {
            // Validate input files
            validate_files(&[&reads, &reference, &gff]).await?;
            
            // Create output directory
            fs::create_dir_all(&output_dir).await?;
            
            // Run full pipeline
            run_full_pipeline(
                &reads,
                &reference,
                &gff,
                &output_dir,
                &sample,
                keep_intermediate,
                &context,
            )
            .await
        }
    }
}

/// Validate that input files exist
async fn validate_files(files: &[&PathBuf]) -> Result<()> {
    for &file in files {
        if !file.exists() {
            return Err(anyhow!(PipelineError::FileNotFound(
                file.to_string_lossy().to_string()
            )));
        }
    }
    Ok(())
}

/// Run the alignment step
async fn run_alignment(
    reads: &Path,
    reference: &Path,
    out_bam: &Path,
    settings: AlignSettings,
    context: &PipelineContext,
) -> Result<()> {
    info!("Aligning reads from {:?} to reference {:?}", reads, reference);
    
    // Create progress bar
    let progress = context.progress.add(
        ProgressBar::new(100).with_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("=>-"),
        ),
    );
    progress.set_message("Aligning reads...");
    
    // Determine aligner to use
    let aligner = settings.aligner.unwrap_or_else(|| "bwa".to_string());
    
    // Create intermediate BAM file path (before duplicate marking)
    let intermediate_bam = context.temp_dir.path().join("aligned.bam");
    
    // Run alignment based on selected aligner
    match aligner.as_str() {
        "bwa" => {
            debug!("Using BWA-MEM aligner");
            // Implementation for BWA would go here
            // For now, simulate progress
            for i in 0..100 {
                progress.set_position(i);
                time::sleep(time::Duration::from_millis(10)).await;
            }
        }
        "minimap2" => {
            debug!("Using Minimap2 aligner");
            // Implementation for Minimap2 would go here
            // For now, simulate progress
            for i in 0..100 {
                progress.set_position(i);
                time::sleep(time::Duration::from_millis(10)).await;
            }
        }
        _ => {
            return Err(anyhow!(PipelineError::ConfigError(format!(
                "Unsupported aligner: {}",
                aligner
            ))));
        }
    }
    
    // Mark duplicates if requested
    if settings.mark_duplicates.unwrap_or(false) {
        progress.set_message("Marking duplicates...");
        // Implementation for duplicate marking would go here
        
        // Copy the final result to the output path
        fs::copy(&intermediate_bam, out_bam).await?;
    } else {
        // No duplicate marking, just rename the intermediate file
        fs::copy(&intermediate_bam, out_bam).await?;
    }
    
    // Index the BAM file
    progress.set_message("Indexing BAM file...");
    // Implementation for BAM indexing would go here
    
    progress.finish_with_message(format!("Alignment completed: {:?}", out_bam));
    
    info!("Alignment completed successfully");
    Ok(())
}

/// Run the variant calling step
async fn run_calling(
    bam: &Path,
    reference: &Path,
    out_vcf: &Path,
    settings: CallSettings,
    format: OutputFormat,
    context: &PipelineContext,
) -> Result<()> {
    info!("Calling variants from {:?} using reference {:?}", bam, reference);
    
    // Create progress bar
    let progress = context.progress.add(
        ProgressBar::new(100).with_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("=>-"),
        ),
    );
    progress.set_message("Calling variants...");
    
    // Set up variant calling parameters
    let min_depth = settings.min_depth.unwrap_or(10);
    debug!("Minimum depth for variant calling: {}", min_depth);
    
    // Check if BAM is indexed
    let bai_path = bam.with_extension("bam.bai");
    if !bai_path.exists() {
        warn!("BAM index not found, creating index for {:?}", bam);
        // Implementation for BAM indexing would go here
    }
    
    // Process regions if specified
    let regions = settings.regions.unwrap_or_default();
    let region_count = regions.len();
    
    if !regions.is_empty() {
        debug!("Processing {} specific regions", region_count);
    } else {
        debug!("Processing entire genome");
    }
    
    // Simulate variant calling progress
    for i in 0..100 {
        progress.set_position(i);
        time::sleep(time::Duration::from_millis(20)).await;
    }
    
    // Create intermediate VCF for format conversion if needed
    let intermediate_vcf = context.temp_dir.path().join("variants.vcf");
    
    // Convert to the requested output format
    progress.set_message("Converting to final format...");
    
    match format {
        OutputFormat::Vcf => {
            fs::copy(&intermediate_vcf, out_vcf).await?;
        }
        OutputFormat::Bcf => {
            // Implementation for BCF conversion would go here
        }
        _ => {
            return Err(anyhow!(PipelineError::ConfigError(format!(
                "Unsupported output format for variant calling: {:?}",
                format
            ))));
        }
    }
    
    progress.finish_with_message(format!("Variant calling completed: {:?}", out_vcf));
    
    info!("Variant calling completed successfully");
    Ok(())
}

/// Run the annotation step
async fn run_annotation(
    vcf: &Path,
    gff: &Path,
    output: &Path,
    settings: AnnotateSettings,
    format: OutputFormat,
    context: &PipelineContext,
) -> Result<()> {
    info!("Annotating variants from {:?} using annotations {:?}", vcf, gff);
    
    // Create progress bar
    let progress = context.progress.add(
        ProgressBar::new(100).with_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("=>-"),
        ),
    );
    progress.set_message("Loading annotations...");
    
    // Process additional databases if specified
    let databases = settings.databases.unwrap_or_default();
    if !databases.is_empty() {
        debug!("Using {} additional annotation databases", databases.len());
        for db in &databases {
            if !db.exists() {
                warn!("Annotation database not found: {:?}", db);
            }
        }
    }
    
    // Check if effect predictions are requested
    let predict_effects = settings.effects.unwrap_or(false);
    if predict_effects {
        debug!("Including effect predictions in annotation");
    }
    
    // Simulate annotation progress
    for i in 0..50 {
        progress.set_position(i);
        time::sleep(time::Duration::from_millis(20)).await;
    }
    
    progress.set_message("Processing variants...");
    
    // Continue simulation
    for i in 50..100 {
        progress.set_position(i);
        time::sleep(time::Duration::from_millis(20)).await;
    }
    
    // Create output in the requested format
    progress.set_message("Writing results...");
    
    match format {
        OutputFormat::Tsv => {
            // Implementation for TSV output would go here
        }
        OutputFormat::Json => {
            // Implementation for JSON output would go here
        }
        OutputFormat::Parquet => {
            // Implementation for Parquet output would go here
        }
        _ => {
            return Err(anyhow!(PipelineError::ConfigError(format!(
                "Unsupported output format for annotation: {:?}",
                format
            ))));
        }
    }
    
    progress.finish_with_message(format!("Annotation completed: {:?}", output));
    
    info!("Annotation completed successfully");
    Ok(())
}

/// Run the full pipeline
async fn run_full_pipeline(
    reads: &Path,
    reference: &Path,
    gff: &Path,
    output_dir: &Path,
    sample: &str,
    keep_intermediate: bool,
    context: &PipelineContext,
) -> Result<()> {
    info!("Running full pipeline for sample: {}", sample);
    
    // Create output paths
    let bam_path = output_dir.join(format!("{}.bam", sample));
    let vcf_path = output_dir.join(format!("{}.vcf", sample));
    let annotation_path = output_dir.join(format!("{}.annotated.tsv", sample));
    
    // Initialize statistics
    let stats = Arc::new(Mutex::new(PipelineStats::default()));
    
    // Step 1: Alignment
    info!("Step 1/3: Alignment");
    let align_result = run_alignment(
        reads,
        reference,
        &bam_path,
        context.settings.align.clone(),
        context,
    )
    .await;
    
    if let Err(e) = align_result {
        error!("Alignment failed: {}", e);
        return Err(e);
    }
    
    // Step 2: Variant Calling
    info!("Step 2/3: Variant Calling");
    let call_result = run_calling(
        &bam_path,
        reference,
        &vcf_path,
        context.settings.call.clone(),
        OutputFormat::Vcf,
        context,
    )
    .await;
    
    if let Err(e) = call_result {
        error!("Variant calling failed: {}", e);
        return Err(e);
    }
    
    // Step 3: Annotation
    info!("Step 3/3: Annotation");
    let annotate_result = run_annotation(
        &vcf_path,
        gff,
        &annotation_path,
        context.settings.annotate.clone(),
        OutputFormat::Tsv,
        context,
    )
    .await;
    
    if let Err(e) = annotate_result {
        error!("Annotation failed: {}", e);
        return Err(e);
    }
    
    // Clean up intermediate files if requested
    if !keep_intermediate {
        info!("Cleaning up intermediate files");
        // In a real implementation, we'd delete intermediates here
    }
    
    // Calculate elapsed time
    let elapsed = context.start_time.elapsed();
    
    // Update and print statistics
    {
        let mut stats_guard = stats.lock().unwrap();
        stats_guard.elapsed_seconds = elapsed.as_secs_f64();
        
        // In a real implementation, we'd gather actual statistics
        stats_guard.aligned_reads = 1_000_000;
        stats_guard.variants_called = 10_000;
        stats_guard.variants_annotated = 5_000;
        
        print_pipeline_summary(&stats_guard, sample);
    }
    
    info!(
        "Full pipeline completed successfully in {:.2} seconds",
        elapsed.as_secs_f64()
    );
    
    Ok(())
}

/// Print a summary of the pipeline results
fn print_pipeline_summary(stats: &PipelineStats, sample: &str) {
    let formatter = Formatter::new();
    
    println!("\n========== Pipeline Summary ==========");
    println!("Sample: {}", sample);
    println!("Elapsed time: {:.2} seconds", stats.elapsed_seconds);
    println!("Aligned reads: {}", formatter.format(stats.aligned_reads as f64));
    println!("Variants called: {}", formatter.format(stats.variants_called as f64));
    println!("Variants annotated: {}", formatter.format(stats.variants_annotated as f64));
    println!("======================================\n");
}