#!/usr/bin/env python3
import json
import random
import argparse
import os

def generate_random_sequence(length: int) -> str:
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))

def generate_kmers_from_sequence(sequence: str, k: int):
    if len(sequence) < k:
        return []
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def generate_transcript_sequences(num_transcripts: int, min_length: int = 500, max_length: int = 3000):
    transcripts = {}
    for i in range(num_transcripts):
        transcript_id = f"transcript_{i:04d}"
        length = random.randint(min_length, max_length)
        sequence = generate_random_sequence(length)
        transcripts[transcript_id] = sequence
    return transcripts

def generate_kmer_index(transcripts, k: int = 31):
    kmer_to_transcripts = {}
    
    for transcript_id, sequence in transcripts.items():
        kmers = generate_kmers_from_sequence(sequence, k)
        for kmer in kmers:
            if kmer not in kmer_to_transcripts:
                kmer_to_transcripts[kmer] = []
            kmer_to_transcripts[kmer].append(transcript_id)
    
    kmer_index = []
    for kmer, transcript_list in kmer_to_transcripts.items():
        kmer_index.append({
            "kmer": kmer,
            "transcripts": list(set(transcript_list)),
            "transcript_positions": None
        })
    
    return kmer_index

def introduce_errors(sequence: str, error_rate: float) -> str:
    bases = ['A', 'T', 'G', 'C']
    seq_list = list(sequence)
    
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            original_base = seq_list[i]
            new_base = random.choice([b for b in bases if b != original_base])
            seq_list[i] = new_base
    
    return ''.join(seq_list)

def write_fastq(reads, filename: str):
    with open(filename, 'w') as f:
        for i, read in enumerate(reads):
            read_id = f"read_{i:08d}"
            quality = 'I' * len(read)
            
            f.write(f"@{read_id}\n")
            f.write(f"{read}\n")
            f.write(f"+\n")
            f.write(f"{quality}\n")

def generate_expression_profile(num_transcripts: int):
    import math
    
    expression = {}
    for i in range(num_transcripts):
        transcript_id = f"transcript_{i:04d}"
        log_expr = random.normalvariate(2, 1.5)
        expression[transcript_id] = math.exp(log_expr)
    
    return expression

def generate_weighted_reads(transcripts, expression_profile, total_reads: int, read_length: int = 100, error_rate: float = 0.01):
    reads = []
    
    total_expression = sum(expression_profile.values())
    normalized_expression = {
        tid: expr / total_expression 
        for tid, expr in expression_profile.items()
    }
    
    transcript_reads = {}
    for transcript_id, expr_fraction in normalized_expression.items():
        transcript_reads[transcript_id] = int(total_reads * expr_fraction)
    
    for transcript_id, num_reads in transcript_reads.items():
        if transcript_id not in transcripts:
            continue
            
        transcript_seq = transcripts[transcript_id]
        if len(transcript_seq) < read_length:
            continue
        
        for _ in range(num_reads):
            start_pos = random.randint(0, len(transcript_seq) - read_length)
            read_seq = transcript_seq[start_pos:start_pos + read_length]
            
            if error_rate > 0:
                read_seq = introduce_errors(read_seq, error_rate)
            
            reads.append(read_seq)
    
    random.shuffle(reads)
    return reads

def main():
    parser = argparse.ArgumentParser(description='Generate sample data for pseudo-alignment')
    parser.add_argument('--num-transcripts', type=int, default=500, help='Number of transcripts')
    parser.add_argument('--num-reads', type=int, default=50000, help='Number of reads')
    parser.add_argument('--read-length', type=int, default=100, help='Read length')
    parser.add_argument('--kmer-length', type=int, default=31, help='K-mer length')
    parser.add_argument('--error-rate', type=float, default=0.02, help='Sequencing error rate')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    random.seed(args.seed)
    os.makedirs('data', exist_ok=True)
    
    print(f"Generating {args.num_transcripts} transcripts and {args.num_reads} reads...")
    
    transcripts = generate_transcript_sequences(args.num_transcripts)
    expression_profile = generate_expression_profile(args.num_transcripts)
    kmer_index = generate_kmer_index(transcripts, args.kmer_length)
    
    with open('data/kmer_index.json', 'w') as f:
        json.dump(kmer_index, f, indent=2)
    
    reads = generate_weighted_reads(
        transcripts, 
        expression_profile, 
        args.num_reads,
        args.read_length,
        args.error_rate
    )
    
    write_fastq(reads, 'data/reads.fastq')
    
    with open('data/transcripts.fasta', 'w') as f:
        for transcript_id, sequence in transcripts.items():
            f.write(f">{transcript_id}\n{sequence}\n")
    
    with open('data/true_expression.tsv', 'w') as f:
        f.write("transcript_id\ttrue_expression\n")
        for transcript_id, expr in sorted(expression_profile.items()):
            f.write(f"{transcript_id}\t{expr:.6f}\n")
    
    print("✓ Sample dataset generated successfully!")
    print(f"  - {len(kmer_index)} k-mers in index")
    print(f"  - {len(reads)} reads generated")
    print("  - Files: data/kmer_index.json, data/reads.fastq")

if __name__ == "__main__":
    main()
