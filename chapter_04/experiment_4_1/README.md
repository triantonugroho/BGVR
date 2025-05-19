## 4.1. Introduction to Functional Genomics Data Structures

### experiment_41

In this simplified example, we demonstrate a workflow for reading a synthetic FASTQ file to build a Position Weight Matrix (PWM) and then a simple Markov Random Field (MRF) using Rust. A stand alone Nextflow serves as our workflow engine, orchestrating the steps and ensuring reproducibility. Rust’s safety and performance make it an excellent choice for bioinformatics tasks, while Nextflow provides a robust framework for pipeline development and execution across diverse compute environments.

The Nextflow pipeline reads a synthetic FASTQ file (synthetic_reads.fastq) and sends it to a single process that compiles a Rust program via Cargo, then runs the resulting binary. Within this Rust program, each FASTQ record is parsed to retrieve its sequence data, and these sequences—assumed to be of the same length for simplicity—are used to build a Position Weight Matrix (PWM) by counting nucleotide frequencies per position and normalizing them to probabilities. Next, a simple first-order Markov Random Field (MRF) model is constructed by tallying each pair of adjacent nucleotides (e.g., A→C) across all sequences and converting these counts to transition probabilities. Finally, the results are written to two separate files—pwm_results.txt for the position-wise probabilities and mrf_results.txt for the transition probabilities—thereby completing the workflow.

#### Project Structure:

```plaintext
experiment_41/
├── Cargo.toml                        # Rust project configuration and dependencies
└── src/
    ├── main.rs                       # Main Rust script containing program logic
    ├── main.nf                       # Nextflow workflow script
    ├── synthetic_reads.fastq.rar     # Compressed synthetic reads FASTQ file
    ├── mrf_results.txt               # MRF results output file
    └── pwm_results.txt               # PWM results output file
```

#### How to run:

run main.rs in powershell:

```powershell
cargo run --release "C:\Users\trian\BGVR\chapter_04\experiment_41\src\synthetic_reads.fastq" pwm_results.txt mrf_results.txt 
```

run main.nf in WSL:

```wsl
nextflow run main.nf --synthetic_fastq 'synthetic_reads.fastq'
```

#### [dependencies]

```toml
bio = "2.0.3"
```

#### Explanation of the Output

##### 1. PWM Results (pwm_results.txt)
The Position Weight Matrix (PWM) provides the probabilities of observing each nucleotide (A, C, G, T) at specific positions within the sequence.

* Each row corresponds to a specific position in the sequence.
* The values represent the probability of finding a particular nucleotide at that position.
* For example:
  * Position 0: 𝐴 = 0.977, 𝐶 = 0.017, 𝐺 = 0.004, 𝑇 = 0.002
    → This indicates that at position 0, nucleotide A has a 97.7% chance of appearing.
  * Position 16: 𝐴 = 0.001, 𝐶 = 0.002, 𝐺 = 0.002, 𝑇 = 0.995
    → At position 16, nucleotide T is highly dominant with a 99.5% probability.

###### Interpretation:
* High probabilities for specific nucleotides at certain positions suggest sequence conservation or motif presence.
* Positions with evenly distributed probabilities suggest low specificity or variability at that site.

##### 2. MRF Results (mrf_results.txt)
The Markov Random Field (MRF) results represent the transition probabilities between nucleotides.

* Each row shows the probability of transitioning from one nucleotide to another.
* For example:
  * 𝐴 → 𝐴 = 0.3013 → Probability of A being followed by A is 30.13%
  * 𝐺 → 𝑇 = 0.2411 → Probability of G being followed by T is 24.11%

###### Interpretation:
* Higher transition probabilities between certain pairs suggest sequence patterns or dependencies.
* Lower transition probabilities imply that such transitions are rare in the sequence.

#### Conclusion
* The PWM and MRF outputs provide complementary insights into the sequence:
  * PWM shows the positional preferences of nucleotides, which helps identify conserved motifs.
  * MRF shows the transition patterns between nucleotides, revealing underlying sequence structure or biases.
* These results can be used to:
  * Identify binding motifs in DNA sequences.
  * Predict sequence patterns based on transition probabilities.
  * Improve sequence alignment or motif search algorithms.
