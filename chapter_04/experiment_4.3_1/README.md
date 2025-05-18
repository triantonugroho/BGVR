## 4.3. Motif Discovery and Regulatory Element Identification

### experiment_4.3_1

The following Rust snippet demonstrates a minimal hidden Markov model (HMM) for a motif-finding scenario, treating “motif” and “non-motif” as distinct states in a DNA sequence. Each state has specified transition probabilities (how likely it is to remain in the current state or switch to the other) and emission probabilities (how likely each nucleotide is under that state). In a real setup, multiple states might represent different motif positions or more complex background distributions. This toy example focuses on capturing whether a genomic position is likely within a motif or background.

The MotifHMM struct holds transition and emission tables in hash maps. The core logic is in the viterbi method, which calculates a most-likely path of states for an input sequence. It initializes dynamic programming arrays for each position in the sequence, storing the highest log probability of arriving in each state and a backpointer to recover the best path. Iterating forward, it updates log probabilities using the formula

\text{score} = \text{prev\_score} + \ln P(\text{transition}) + \ln P(\text{emission}),

selecting whichever predecessor state yields the greatest value. After processing all positions, the algorithm identifies the best final state and follows backpointers to reconstruct the most-likely motif vs. non-motif labeling of the entire DNA sequence.

#### Project Structure:

```plaintext
experiment_43_1/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_4.3_1"
version = "0.1.0"
edition = "2024"

[dependencies]

```

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```
(run main.rs and save the output in output.txt)


#### Explanation of the Output
The output of the Hidden Markov Model (HMM) motif-finding program shows the input DNA sequence and the most likely state path computed by the Viterbi algorithm. Let's go through the details step-by-step:

##### 1. Input DNA Sequence

```rust
Sequence:     ACGGAATACACGG
```

* The input sequence is ACGGAATACACGG, which consists of 13 nucleotides.
* The Viterbi algorithm is tasked with identifying whether each position in the sequence is in a Motif or NonMotif state.

#### 2. Most-Likely Path of States

```rust
Most-likely:  [NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif, NonMotif]
```

* The output shows that the most likely state for every position in the sequence is NonMotif.
* This means the algorithm predicts that none of the positions in the sequence are part of a motif.

##### Why Did It Predict All NonMotif?

###### 1. Transition Probabilities
The transition probabilities are set as:

* 𝑃(Motif → Motif) = 0.8
* 𝑃(Motif → NonMotif) = 0.2
* 𝑃(NonMotif → Motif) = 0.05
* 𝑃(NonMotif → NonMotif) = 0.95

This means the HMM is biased toward staying in the NonMotif state since 𝑃(NonMotif → NonMotif) = 0.95 is very high.

##### 2. Emission Probabilities
The emission probabilities are set as:

* NonMotif State: Equal probability for all nucleotides (A, C, G, T) = 0.25
* Motif State: Higher probability for A = 0.4 (other nucleotides = 0.2)
Since the emission probabilities for NonMotif are fairly balanced, the model is unlikely to switch to the Motif state unless the sequence contains a clear motif pattern dominated by A.

##### 3. Initial Probability and Bias
The initial state probabilities are equally split between Motif and NonMotif (0.5 each).
However, the model quickly converges to the NonMotif state because the transition and emission biases favor NonMotif.

#### Conclusion
* The HMM and Viterbi algorithm executed successfully.
* The algorithm predicted that all positions are in the NonMotif state due to strong bias in the transition and emission probabilities.
* To improve motif detection:
  * Lower the probability of staying in the NonMotif state.
  * Adjust the emission probabilities for the Motif state to make motifs more distinguishable.
  * Test on sequences known to contain motifs to adjust parameters.

