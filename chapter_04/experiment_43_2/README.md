## 4.3. Motif Discovery and Regulatory Element Identification

### experiment_43_2

The code defines a MotifModel with a PWM representing the current motif guess. In the em_iteration function, each HPC node performs a local E-step by evaluating how well each DNA sequence aligns with the motif, accumulating partial counts in a local array. These partial counts are then merged (e.g., via MPI) into a global count matrix, on which the M-step is applied to recalculate PWM entries. Although real MEME has more complex components (like handling multiple motifs, background models, or variable motif lengths), this toy version shows the core dynamic: partial responsibilities are computed locally, then combined to update the motif model. Rust’s concurrency and HPC-oriented structure can scale this approach to large sets of DNA reads.

This code simulates a single iteration of an Expectation-Maximization (EM)-style motif-finding approach in a high-performance computing setting. The MotifModel struct holds a position weight matrix (PWM), initialized randomly in new_random. The em_iteration function (the E-step) processes a subset of sequences by sliding the PWM along each read, selecting the highest-probability alignment, and tallying partial counts in a local partial_counts array. In a real HPC pipeline, multiple nodes would each produce such partial counts, which would then be merged (e.g., by summing across nodes). Finally, model.update_from_counts (the M-step) re-estimates PWM entries to reflect this iteration’s aggregated partial data. Over multiple iterations—and across multiple HPC nodes—this procedure converges on a refined motif model that best explains the data segments distributed across the cluster.

#### Files contents:
* experiment_43_2/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_43_2/src/
  * main.rs (rust script)
  * output.txt (output file)

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```
(run main.nf and save the output in output.txt)

#### [dependencies]

```toml
rand = "0.9.0"
```

#### Explanation of the Output
The program is implementing an Expectation-Maximization (EM) algorithm for motif finding using a Position Weight Matrix (PWM). Let's walk through the output:

###### Initial Motif Model

```rust
Initial motif model: [[0.30213606463566717, 0.3665743064611263, 0.3257436654557016, 0.005545963447505091], 
                      [0.11302343853582049, 0.21140046747006883, 0.35164299943680016, 0.3239330945573104], 
                      [0.21824434961237274, 0.009446313565084886, 0.5070345041590435, 0.26527483266349894], 
                      [0.3209274488558195, 0.27610623878000706, 0.27159756895679543, 0.13136874340737817]]
```

##### What this means:
* The motif model is initialized randomly using new_random().
* Each row corresponds to a position in the motif.
* Each value in a row is the probability of observing one of the four nucleotides (A, C, G, T) at that position.
* For example:
  * At position 1, the probability of observing A, C, G, and T is approximately [0.302, 0.367, 0.326, 0.006].

##### Updated Motif Model

```rust
Updated motif model: [[0.0, 0.3333333333333333, 0.6666666666666666, 0.0], 
                      [0.0, 0.0, 0.6666666666666666, 0.3333333333333333], 
                      [0.0, 0.0, 0.3333333333333333, 0.6666666666666666], 
                      [0.6666666666666666, 0.0, 0.3333333333333333, 0.0]]
```

##### What happened:

###### 1. After the E-step and M-step:

* The model updates the PWM based on the aligned sequences.
* The algorithm identifies the most likely motif positions and updates the counts.
* These counts are normalized into probabilities in update_from_counts().

##### The model converges to a cleaner pattern:

* First position: G dominates with 66% probability, C follows with 33%.
* Second position: G still dominates.
* Third position: T becomes dominant.
* Fourth position: A becomes the most probable nucleotide.

##### Interpretation
* The EM algorithm is working as expected:
* In the E-step: The model identifies the best motif alignment in each sequence.
* In the M-step: The model adjusts the PWM using partial counts.
* The PWM has converged to reflect the motif patterns observed in the sequences.
* The updated PWM is sharper (more polarized probabilities), suggesting the model has successfully learned a motif pattern from the data.

#### Conclusion
* The EM algorithm executed successfully.
* The PWM learned a clear motif pattern from the input sequences.
* The motif is represented by high probabilities in certain nucleotide positions, indicating that the model has focused on consistent patterns.
* The model is ready for further refinement or testing on larger datasets.
