use anyhow::{Context, Result};
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::fs::{File, create_dir_all};
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::sync::Arc;
use clap::Parser;

#[derive(Serialize, Deserialize, Debug)]
struct CorrectionPatch {
    contig: String,
    position: u64,
    ref_base: u8,
    suggested_base: u8,
    coverage: u64,
}

#[derive(Parser, Debug)]
#[command(name = "polish_patches")]
struct Args {
    #[arg(long = "bam-input", alias = "bam")]
    bam_input: PathBuf,
    #[arg(long)]
    gtf: PathBuf, // Tambahkan ini untuk mendukung argumen --gtf
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = 10_000)]
    chunk_size: usize,
    #[arg(long, default_value = "partial_patches")]
    partial_outdir: PathBuf,
    #[arg(long, default_value = "merged_corrections.json")]
    merged_output: PathBuf,
}

#[derive(Serialize, Deserialize, Debug)]
struct PartialPatches {
    patches: Vec<CorrectionPatch>,
}

fn decode_base(encoded: u8) -> u8 {
    match encoded {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}

fn read_chunk(bam_reader: &mut bam::Reader, chunk_size: usize) -> Result<Vec<bam::Record>> {
    let mut chunk = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        if let Some(result) = bam_reader.records().next() {
            chunk.push(result?);
        } else {
            break;
        }
    }
    Ok(chunk)
}

fn prepare_contig_map(header: &bam::HeaderView) -> Arc<HashMap<i32, String>> {
    let contig_map: HashMap<i32, String> = (0..header.target_count() as i32)
        .map(|tid| (tid, String::from_utf8_lossy(header.tid2name(tid as u32)).to_string()))
        .collect();
    Arc::new(contig_map)
}

fn collect_correction_patches(records: &[bam::Record], contig_map: Arc<HashMap<i32, String>>) -> Vec<CorrectionPatch> {
    let mut grouped_records: HashMap<i32, Vec<bam::Record>> = HashMap::new();
    
    for record in records {
        grouped_records.entry(record.tid()).or_default().push(record.clone());
    }

    grouped_records.into_par_iter().flat_map(|(tid, recs)| {
        let contig_name = contig_map.get(&tid).unwrap_or(&"unknown".to_string()).clone();
        let mut mismatch_map: HashMap<u64, (u8, u8, u64)> = HashMap::new();

        for r in recs {
            let start_pos = r.pos() as u64;
            let seq = r.seq();
            for (i, base_enc) in seq.as_bytes().iter().enumerate() {
                let global_pos = start_pos + i as u64;
                let read_base = decode_base(*base_enc);
                if read_base != b'N' {
                    if read_base != b'A' {
                        let entry = mismatch_map.entry(global_pos).or_insert((b'A', read_base, 0));
                        entry.2 += 1;
                    }
                }
            }
        }

        mismatch_map.into_iter().map(|(pos, (ref_base, suggested, coverage))| {
            CorrectionPatch {
                contig: contig_name.clone(),
                position: pos,
                ref_base,
                suggested_base: suggested,
                coverage,
            }
        }).collect::<Vec<_>>()
    }).collect()
}

fn merge_patch_vectors(mut acc: Vec<CorrectionPatch>, mut new_patches: Vec<CorrectionPatch>) -> Vec<CorrectionPatch> {
    acc.append(&mut new_patches);
    acc
}

fn main() -> Result<()> {
    let args = Args::parse();
    create_dir_all(&args.partial_outdir)?;

    let mut bam_reader = bam::Reader::from_path(&args.bam_input)?;
    let contig_map = prepare_contig_map(&bam_reader.header());

    let mut chunk_index = 0;
    loop {
        let chunk = read_chunk(&mut bam_reader, args.chunk_size)?;
        if chunk.is_empty() {
            break;
        }

        let partial_patches = collect_correction_patches(&chunk, Arc::clone(&contig_map));
        let container = PartialPatches { patches: partial_patches };

        let outfile_path = args.partial_outdir.join(format!("partial_patches_chunk_{}.json", chunk_index));
        let out_file = File::create(&outfile_path)?;
        serde_json::to_writer(BufWriter::new(out_file), &container)?;

        println!("Processed chunk {} with {} records. Partial patches stored at {:?}.",
            chunk_index, chunk.len(), outfile_path);
        chunk_index += 1;
    }

    let dir_entries = std::fs::read_dir(&args.partial_outdir)?;
    let mut merged_patches = Vec::new();
    
    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |p| p.to_string_lossy().starts_with("partial_patches_chunk_")) {
            let file = File::open(&path)?;
            let partial: PartialPatches = serde_json::from_reader(BufReader::new(file))?;
            merged_patches = merge_patch_vectors(merged_patches, partial.patches);
        }
    }

    let merged_file = File::create(&args.merged_output)?;
    serde_json::to_writer(BufWriter::new(merged_file), &merged_patches)?;

    println!("Merged {} total correction patches into {:?}.", merged_patches.len(), args.merged_output);
    Ok(())
}
