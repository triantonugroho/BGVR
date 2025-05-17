## 4.3. Motif Discovery and Regulatory Element Identification

### experiment_43_3

This Rust code shows how one might implement a parallelized Gibbs sampling routine for motif discovery. Each parallel chain maintains a separate “guessed” motif configuration, sampling new motif positions from a conditional probability distribution based on the other sequences’ motif assignments. In a real HPC environment, multiple processes or threads each run a Gibbs chain independently, possibly seeded with different random initial states, and after a specified number of iterations, the partial results are consolidated (e.g., by picking the chain with the highest posterior or by averaging motif parameters).

In the sample code, the GibbsSampler struct holds a motif model (like a PWM) and the sequences. The run_one_iteration method picks a random sequence, “unassigns” its motif position, then samples a new position for that sequence’s motif from a distribution proportional to the motif’s likelihood. Real Bayesian approaches would also maintain posterior distributions for the motif parameters, but for brevity this snippet only tracks the motif location in each sequence. The run_parallel_chains function demonstrates launching multiple threads (rayon::scope) where each chain runs independently for a fixed number of iterations. After the chains complete, one might compare or merge the results. This design captures the central idea of parallelizing Gibbs sampling: each chain is a Markov chain that explores motif placements, and concurrency reduces overall runtime for large data sets.

#### Project Structure:

```plaintext
experiment_43_3/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Output file
``

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```
(run main.rs and save the output in output.txt)

#### [dependencies]

```toml
rand = "0.9.0"
rayon = "1.10.0"
```

#### Explanation of Output
The main.rs code implements a Gibbs sampling algorithm for motif discovery in DNA sequences. Here's a step-by-step breakdown of what the output represents:

##### 1. Input Data

* sequences: Three DNA sequences are provided as input:
```rust
let sequences = vec![
    b"ACGATGATGAC".to_vec(),   // Length = 11
    b"TTTTAAAACCCCGG".to_vec(), // Length = 14
    b"AAAATGATGAAAA".to_vec(),  // Length = 13
];
```

* k = 5: The motif length is 5 (the length of the subsequence that the algorithm tries to align across sequences).
* num_chains = 3: Three independent Gibbs sampling chains are executed in parallel.
* iterations = 10: Each chain performs 10 iterations of motif refinement.

##### 2. Initial Setup

* The algorithm initializes the starting position of motifs randomly in each sequence.
* For each chain, a GibbsSampler instance is created with random starting positions for the motifs.

##### 3. Sampling Process

* In each Gibbs sampling iteration:
  * A random sequence is selected.
  * The motif position for that sequence is "unassigned."
  * A new position is sampled based on the computed likelihood of each possible k-mer (subsequence of length k).
  * The likelihood function is simplistic — it assigns higher scores to k-mers containing more A or a characters.
* After each iteration, the motif position for the selected sequence is updated based on the sampling outcome.

##### 4. Parallel Execution

* The chains are run in parallel using rayon::scope, ensuring independent execution of each Gibbs sampling chain.
* The results from each chain are collected using an Arc<Mutex<>> structure, allowing safe concurrent access.

##### 5. Final Results

After 10 iterations of refinement in each chain, the final motif positions are printed:

###### Output:

```rust
Sampler 0 final motif positions: [3, 4, 4]
Sampler 1 final motif positions: [5, 5, 0]
Sampler 2 final motif positions: [3, 5, 5]
```

#### Interpretation of Results
Each array corresponds to the final motif start positions in the three sequences, output by each of the three chains:

* Sampler 0:

  * Sequence 1 → Starting at position 3
  * Sequence 2 → Starting at position 4
  * Sequence 3 → Starting at position 4

* Sampler 1:

  * Sequence 1 → Starting at position 5
  * Sequence 2 → Starting at position 5
  * Sequence 3 → Starting at position 0

* Sampler 2:

  * Sequence 1 → Starting at position 3
  * Sequence 2 → Starting at position 5
  * Sequence 3 → Starting at position 5

#### Conclusion
* The Gibbs sampler successfully:
  * Initialized motifs randomly.
  * Refined motif positions over 10 iterations based on likelihood.
  * Executed multiple chains in parallel without conflicts.
  * Final motif positions reflect the convergence of the algorithm toward positions with higher scores (more A characters).
