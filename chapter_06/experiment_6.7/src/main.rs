// Replace log crate macros with println!
use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::bcf::{Reader as BcfReader, Read as BcfRead};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

/// A structure to store basic variant information, including coverage and optional annotation.
#[derive(Debug)]
struct VariantInfo {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    coverage: u32,
    annotation: Option<String>,
}

/// Command-line arguments for the integrative analysis tool.
#[derive(Parser, Debug)]
#[command(name = "rust_integrative_tool", about = "A robust Rust-based tool for coverage and variant annotation")]
struct Cli {
    /// Path to the indexed BAM file
    #[arg(long, value_name = "BAM_FILE")]
    bam: String,

    /// Path to a BCF or VCF file containing variant calls (must be indexed)
    #[arg(long, value_name = "BCF_FILE")]
    bcf: String,

    /// Path to a GFF annotation file (converted to a simple HashMap for demonstration)
    #[arg(long, value_name = "GFF_ANNOTATION")]
    gff: String,
}

fn load_gff_annotations(gff_path: &str) -> Result<HashMap<(String, u64, u64), String>> {
    let mut annotations = HashMap::new();
    let content = std::fs::read_to_string(gff_path)
        .with_context(|| format!("Failed to read GFF annotation file: {gff_path}"))?;
    
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let chrom = cols[0].to_string();
        let start: u64 = cols[3].parse().with_context(|| "Parsing start coordinate")?;
        let end: u64 = cols[4].parse().with_context(|| "Parsing end coordinate")?;
        let annotation = cols[8].to_string();
        
        annotations.insert((chrom, start, end), annotation);
    }
    Ok(annotations)
}

fn process_variant(
    bam_path: &str,
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    gff_anno: &HashMap<(String, u64, u64), String>,
) -> Result<VariantInfo> {
    let region_string = format!("{}:{}-{}", chrom, pos + 1, pos + 1);
    let mut bam_reader = IndexedReader::from_path(Path::new(bam_path))
        .with_context(|| format!("Failed to open indexed BAM file: {}", bam_path))?;
    bam_reader.fetch(&region_string)
        .with_context(|| format!("Failed to fetch region: {}", region_string))?;

    let mut coverage_count = 0;
    for read_result in bam_reader.records() {
        let _read = read_result.with_context(|| "Error reading a BAM record")?;
        coverage_count += 1;
    }

    let annotation = gff_anno.iter().find_map(|((c, start, end), feat)| {
        if c == chrom && *start <= pos && *end >= pos {
            Some(feat.clone())
        } else {
            None
        }
    });

    Ok(VariantInfo {
        chrom: chrom.to_string(),
        pos,
        ref_allele: ref_allele.to_string(),
        alt_allele: alt_allele.to_string(),
        coverage: coverage_count,
        annotation,
    })
}

fn integrate_alignment_variant_data(
    bam_path: &str,
    bcf_path: &str,
    gff_anno: &HashMap<(String, u64, u64), String>,
) -> Result<Vec<VariantInfo>> {
    let mut bcf_reader = BcfReader::from_path(Path::new(bcf_path))
        .with_context(|| format!("Failed to open BCF file: {bcf_path}"))?;

    let header = bcf_reader.header().clone();

    let mut variants = Vec::new();
    for result_record in bcf_reader.records() {
        let record = result_record.with_context(|| "Error reading a BCF record")?;
        let rid = record
            .rid()
            .ok_or_else(|| anyhow::anyhow!("BCF record has no valid reference ID"))?;
        
        let chrom = String::from_utf8_lossy(header.rid2name(rid)?).to_string();
        let pos = record.pos() as u64;

        let alleles = record.alleles();
        if alleles.len() < 2 {
            continue;
        }
        let ref_allele = String::from_utf8_lossy(alleles[0]).into_owned();
        let alt_allele = String::from_utf8_lossy(alleles[1]).into_owned();

        variants.push((chrom, pos, ref_allele, alt_allele));
    }

    let variant_infos: Vec<_> = variants
        .par_iter()
        .map(|(chrom, pos, ref_allele, alt_allele)| {
            process_variant(
                bam_path,
                chrom,
                *pos,
                ref_allele,
                alt_allele,
                gff_anno,
            )
        })
        .collect();

    let mut results = Vec::with_capacity(variant_infos.len());
    for vinfo in variant_infos {
        match vinfo {
            Ok(info) => results.push(info),
            Err(e) => {
                println!("Error processing variant: {:?}", e);
            }
        }
    }

    Ok(results)
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    println!("Loading GFF annotations from {}", cli.gff);
    let gff_anno = load_gff_annotations(&cli.gff)?;

    println!("Starting integrative analysis on BAM: {}, BCF: {}", cli.bam, cli.bcf);
    let results = integrate_alignment_variant_data(&cli.bam, &cli.bcf, &gff_anno)?;
    println!("Processed {} variants.", results.len());

    for variant in &results {
        println!(
            "Variant: {}:{} ref={} alt={} coverage={} annotation={:?}",
            variant.chrom, variant.pos, variant.ref_allele, variant.alt_allele, variant.coverage, variant.annotation
        );
    }

    Ok(())
}