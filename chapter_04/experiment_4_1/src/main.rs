use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use bio::io::fastq;
use std::collections::HashMap;

fn main() {
    // Expect command line arguments:
    // 1) input FASTQ
    // 2) PWM output file
    // 3) MRF output file
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <input_fastq> <pwm_output> <mrf_output>", args[0]);
        std::process::exit(1);
    }

    let input_fastq = &args[1];
    let pwm_output_path = &args[2];
    let mrf_output_path = &args[3];

    // Read sequences from FASTQ
    let mut seqs = Vec::new();
    let reader = fastq::Reader::from_file(input_fastq).expect("Could not open FASTQ file");
    for record in reader.records() {
        let rec = record.expect("Error reading record");
        seqs.push(rec.seq().to_vec());
    }

    // Convert sequence bytes to String
    let string_seqs: Vec<String> = seqs.iter()
                                       .map(|s| String::from_utf8_lossy(s).into_owned())
                                       .collect();

    // Find the minimum sequence length
    let min_length = string_seqs.iter().map(|s| s.len()).min().unwrap();

    // Ensure all sequences have the same length by truncating longer ones
    let trimmed_seqs: Vec<String> = string_seqs.iter()
        .map(|s| s.chars().take(min_length).collect())  
        .collect();

    println!("All sequences truncated to length: {}", min_length);

    // Build PWM and MRF
    let pwm = build_pwm(&trimmed_seqs);
    let mrf = build_mrf(&trimmed_seqs);

    // Write PWM results
    let pwm_file = File::create(pwm_output_path).expect("Cannot create PWM output file");
    let mut pwm_writer = BufWriter::new(pwm_file);
    writeln!(pwm_writer, "Position Weight Matrix (probabilities)").unwrap();
    for (pos, position_map) in pwm.iter().enumerate() {
        writeln!(
            pwm_writer,
            "Position {}: A={:.3}, C={:.3}, G={:.3}, T={:.3}",
            pos,
            position_map[&'A'],
            position_map[&'C'],
            position_map[&'G'],
            position_map[&'T']
        ).unwrap();
    }

    // Write MRF results
    let mrf_file = File::create(mrf_output_path).expect("Cannot create MRF output file");
    let mut mrf_writer = BufWriter::new(mrf_file);
    writeln!(mrf_writer, "1st-order Markov Random Field (transition probabilities)").unwrap();
    for base1 in &['A', 'C', 'G', 'T'] {
        for base2 in &['A', 'C', 'G', 'T'] {
            let probability = mrf.get(&(*base1, *base2)).copied().unwrap_or(0.0);
            writeln!(mrf_writer, "{}->{}: {:.4}", base1, base2, probability).unwrap();
        }
    }

    println!("PWM and MRF computation completed successfully!");
}

// Function to build PWM
fn build_pwm(seqs: &[String]) -> Vec<HashMap<char, f64>> {
    let seq_len = seqs[0].len();
    let mut pwm = vec![HashMap::new(); seq_len];

    for pos in 0..seq_len {
        let mut counts = HashMap::new();
        for base in &['A', 'C', 'G', 'T'] {
            counts.insert(*base, 0.0); // Initialize counts to avoid missing bases
        }

        for seq in seqs {
            if let Some(base) = seq.chars().nth(pos) {
                if let Some(count) = counts.get_mut(&base) {
                    *count += 1.0;
                }
            }
        }

        // Normalize probabilities
        let total = seqs.len() as f64;
        for base in &['A', 'C', 'G', 'T'] {
            let prob = counts.get(base).copied().unwrap_or(0.0) / total;
            pwm[pos].insert(*base, prob);
        }
    }

    pwm
}


// Function to build MRF
fn build_mrf(seqs: &[String]) -> HashMap<(char, char), f64> {
    let mut transitions = HashMap::new();
    let mut total_counts = HashMap::new();

    // Initialize all base transitions
    for base1 in &['A', 'C', 'G', 'T'] {
        for base2 in &['A', 'C', 'G', 'T'] {
            transitions.insert((*base1, *base2), 0.0);
        }
        total_counts.insert(*base1, 0.0);
    }

    for seq in seqs {
        let chars: Vec<char> = seq.chars().collect();
        for i in 0..(chars.len() - 1) {
            let pair = (chars[i], chars[i + 1]);
            if let Some(count) = transitions.get_mut(&pair) {
                *count += 1.0;
            }
            if let Some(count) = total_counts.get_mut(&chars[i]) {
                *count += 1.0;
            }
        }
    }

    // Normalize probabilities
    for ((base1, base2), count) in transitions.iter_mut() {
        let total = total_counts.get(base1).copied().unwrap_or(1.0);
        *count /= total;
    }

    transitions
}

