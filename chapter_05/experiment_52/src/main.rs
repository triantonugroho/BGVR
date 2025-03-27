use needletail::parser::SequenceRecord;
use needletail::parse_fastx_file;
use std::borrow::Cow;
use std::fs::File;
use std::io::{self, Write};

fn sliding_window_trim(seq: &[u8], quals: &[u8], window_size: usize, quality_threshold: u8) -> (usize, usize) {
    if seq.is_empty() || quals.is_empty() || seq.len() != quals.len() {
        return (0, 0);
    }

    let mut best_start = 0;
    let mut best_end = seq.len();
    let mut current_quality_sum = 0;
    let mut best_quality_sum = 0;

    for i in 0..seq.len() {
        current_quality_sum += quals[i] as usize;
        if i >= window_size {
            current_quality_sum -= quals[i - window_size] as usize;
        }
        if i >= window_size - 1 && current_quality_sum >= (quality_threshold as usize * window_size) {
            best_end = i + 1;
            best_quality_sum = current_quality_sum;
        }
    }

    (best_start, best_end)
}

fn process_chunk(records: &[SequenceRecord], output_prefix: &str, chunk_id: usize, window_size: usize, quality_threshold: u8) -> io::Result<()> {
    let output_filename = format!("{}_chunk_{}.fastq", output_prefix, chunk_id);
    let mut output_file = File::create(output_filename)?;

    for record in records {
        let seq = record.seq().as_ref();
        let quals = record.qual().as_ref();
        let (start, end) = sliding_window_trim(seq, quals, window_size, quality_threshold);

        if start < end {
            writeln!(output_file, "@{}", record.id())?;
            writeln!(output_file, "{}", std::str::from_utf8(&seq[start..end]).unwrap())?;
            writeln!(output_file, "+")?;
            writeln!(output_file, "{}", std::str::from_utf8(&quals[start..end]).unwrap())?;
        }
    }

    Ok(())
}

fn main() {
    let filename = "example.fastq";
    let mut reader = parse_fastx_file(filename).expect("Failed to open file");

    let mut records_buffer: Vec<SequenceRecord> = Vec::new();
    let chunk_size = 100; // Ubah sesuai kebutuhan
    let output_prefix = "output";
    let window_size = 5;
    let quality_threshold = 20;
    let mut chunk_counter = 0;

    while let Some(Ok(record)) = reader.next() {
        records_buffer.push(record);

        if records_buffer.len() >= chunk_size {
            process_chunk(&records_buffer, output_prefix, chunk_counter, window_size, quality_threshold).unwrap();
            records_buffer.clear();
            chunk_counter += 1;
        }
    }

    if !records_buffer.is_empty() {
        process_chunk(&records_buffer, output_prefix, chunk_counter, window_size, quality_threshold).unwrap();
    }
}
