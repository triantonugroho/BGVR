use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, BufReader};
use std::path::PathBuf;
use clap::Parser;

#[derive(Debug, Serialize, Deserialize, Clone)]
struct VariantSite {
    chrom: String,
    position: u64,
    ref_base: char,
    alt_base: char,
    likelihood: f64,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct BaseQuality {
    base: char,
    quality: u8,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct VariantHypothesis {
    chrom: String,
    position: u64,
    ref_base: char,
    alt_base: char,
}

#[derive(Parser, Debug)]
#[command(name = "naive_variant_caller")]
struct Args {
    #[arg(long)]
    pileup_input: PathBuf,
    
    #[arg(long)]
    hypotheses_input: PathBuf,
    
    #[arg(long, default_value_t = 10_000)]
    chunk_size: usize,
    
    #[arg(long, default_value = "partial_variants")]
    output_dir: PathBuf,
    
    #[arg(long, default_value = "merged_variants.json")]
    merged_output: PathBuf,
}

#[derive(Debug, Serialize, Deserialize)]
struct PartialResult {
    sites: Vec<VariantSite>,
}

fn naive_likelihood(reads: &[BaseQuality], ref_base: char, alt_base: char) -> f64 {
    let mut total_ll = 0.0;
    for baseq in reads {
        let p_error = 10f64.powf(-(baseq.quality as f64 / 10.0));
        if baseq.base == alt_base || baseq.base == ref_base {
            total_ll += (1.0 - p_error).ln();
        } else {
            total_ll += p_error.ln();
        }
    }
    total_ll
}

fn compute_variants(
    pileups: &HashMap<(String, u64), Vec<BaseQuality>>,
    hypotheses: &[VariantHypothesis],
) -> Vec<VariantSite> {
    hypotheses
        .par_iter()
        .map(|hyp| {
            let key = (hyp.chrom.clone(), hyp.position);
            let empty_vec = Vec::new();
            let reads = pileups.get(&key).unwrap_or(&empty_vec);

            let ll = naive_likelihood(reads, hyp.ref_base, hyp.alt_base);
            VariantSite {
                chrom: hyp.chrom.clone(),
                position: hyp.position,
                ref_base: hyp.ref_base,
                alt_base: hyp.alt_base,
                likelihood: ll,
            }
        })
        .collect()
}

fn parse_pileup_keys(pileup_raw: HashMap<String, Vec<BaseQuality>>) -> HashMap<(String, u64), Vec<BaseQuality>> {
    let mut pileups = HashMap::new();
    for (key, value) in pileup_raw {
        if let Some((chrom, pos_str)) = key.split_once(':') {
            if let Ok(pos) = pos_str.parse::<u64>() {
                pileups.insert((chrom.to_string(), pos), value);
            }
        }
    }
    pileups
}

fn merge_variant_sites(mut all_sites: Vec<VariantSite>, mut new_sites: Vec<VariantSite>) -> Vec<VariantSite> {
    all_sites.append(&mut new_sites);
    all_sites
}

fn main() -> Result<()> {
    let args = Args::parse();

    create_dir_all(&args.output_dir)?;

    let pileup_file = File::open(&args.pileup_input)?;
    let pileup_raw: HashMap<String, Vec<BaseQuality>> = serde_json::from_reader(BufReader::new(pileup_file))?;
    let pileups = parse_pileup_keys(pileup_raw);

    let hyp_file = File::open(&args.hypotheses_input)?;
    let all_hypotheses: Vec<VariantHypothesis> = serde_json::from_reader(BufReader::new(hyp_file))?;

    let mut start_index = 0;
    let mut chunk_count = 0;
    while start_index < all_hypotheses.len() {
        let end_index = (start_index + args.chunk_size).min(all_hypotheses.len());
        let chunk = &all_hypotheses[start_index..end_index];

        let chunk_sites = compute_variants(&pileups, chunk);
        let partial_res = PartialResult { sites: chunk_sites };
        let chunk_path = args.output_dir.join(format!("partial_variants_chunk_{}.json", chunk_count));
        let out_file = File::create(&chunk_path)?;
        serde_json::to_writer(BufWriter::new(out_file), &partial_res)?;

        println!(
            "Chunk {} processed ({} hypotheses). Partial results saved at {:?}",
            chunk_count, chunk.len(), chunk_path
        );
        start_index = end_index;
        chunk_count += 1;
    }

    let mut merged_variants = Vec::new();
    let dir_entries = std::fs::read_dir(&args.output_dir)?;

    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |p| p.to_string_lossy().starts_with("partial_variants_chunk_")) {
            let file = File::open(&path)?;
            let partial: PartialResult = serde_json::from_reader(BufReader::new(file))?;
            merged_variants = merge_variant_sites(merged_variants, partial.sites);
        }
    }

    let merged_file = File::create(&args.merged_output)?;
    serde_json::to_writer(BufWriter::new(merged_file), &merged_variants)?;

    println!("Successfully merged {} variants into {:?}", merged_variants.len(), args.merged_output);
    Ok(())
}