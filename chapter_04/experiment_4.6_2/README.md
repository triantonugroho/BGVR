## 4.6. Single-Cell Functional Genomics

### experiment_4.6_2

This Rust code provides a high-performance, parallel implementation of sparse matrix–vector multiplication using the Compressed Sparse Row (CSR) format. By leveraging the Rayon library, each row of the CSR matrix is processed concurrently, which makes the code well-suited for large-scale numerical or HPC applications. It demonstrates how Rust’s zero-cost abstractions, compile-time safety checks, and ergonomic concurrency tooling can effectively handle data-intensive workloads.

The core data structure, CsrMatrix, stores non-zero values (values), their corresponding column indices (col_indices), and a set of row pointers (row_ptrs) where row_ptrs[i] marks the start of row i within values and col_indices. In the mul_vector_parallel method, each row is assigned to a different thread through into_par_iter(), and the partial sums (dot products between row segments and the input vector) are computed independently. By accumulating row results in a lock-free way, contention is minimized, and every row’s contribution can be swiftly combined into a final output vector. This approach is highly scalable, allowing larger matrices and vectors to be processed in parallel, and can be further adapted for distributed memory systems or ephemeral container usage in HPC or cloud environments. With Rust’s memory-safety guarantees, large consortia can confidently scale up their single-cell analyses, preserving correctness even in massively parallel environments (Wolf et al. (2018)).

#### Project Structure:

```plaintext
experiment_4.6_2/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Text output file
```

#### Cargo.toml

```toml
[package]
name = "experiment_4.6_2"
version = "0.1.0"
edition = "2024"

[dependencies]
rayon = "1.10.0"
```

#### How to run:

run main.rs in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)


#### Explanation of the Output
My Rust program performs parallel sparse matrix-vector multiplication using the Compressed Sparse Row (CSR) format. The multiplication is done efficiently using Rayon for parallelization and Mutex for thread-safe accumulation.

##### Output Breakdown
The given CSR matrix is:

```rust
    5.0  0.0  0.0  0.0  
    0.0  0.0  2.0  0.0  
    0.0  0.0  0.0  3.0
```
  
* values = [5.0, 2.0, 3.0] → Non-zero elements.
* col_indices = [0, 2, 3] → Corresponding column indices.
* row_ptrs = [0, 1, 2, 3] → Row boundaries in values.

The dense vector is:

```rust
    [2.0, 0.0, 1.0, 1.5]
```

Now, computing CSR * vector row by row:

##### 1. Row 0:

* 5.0 * vec[0] = 5.0 * 2.0 = 10.0
* Result: 10.0

##### 2. Row 1:

* 2.0 * vec[2] = 2.0 * 1.0 = 2.0
* Result: 2.0

##### 3. Row 2:

* 3.0 * vec[3] = 3.0 * 1.5 = 4.5
* Result: 4.5

##### Final output:

```rust
CSR * vector = [10.0, 2.0, 4.5]
```

#### Conclusion
The program correctly computes sparse matrix-vector multiplication using CSR format.

* Parallelization with Rayon improves efficiency, especially for large sparse matrices.
* The Mutex-protected shared result ensures thread safety.
* CSR format saves memory by storing only non-zero values.


