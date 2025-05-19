## 5.1. Introduction to Sequence Analysis

### experiment_51_1

The following code showcases a HyperLogLog (HLL) implementation in Rust, a probabilistic data structure designed for estimating the cardinality (the count of unique elements) of a data set. Unlike traditional methods that store every distinct element, HyperLogLog employs a fixed-size array of small counters, making it extremely memory efficient, particularly for very large data streams. The structure is defined by a parameter p, which determines the number of registers m=2^p. Each incoming item is hashed to a 64-bit value, and the top p bits index into one of these registers. The remaining bits are then examined to find the longest run of leading zeros, and the register is updated with the maximum such run observed over time. By only storing these small counters, HyperLogLog achieves a compact representation of a potentially massive set.

A key feature of HyperLogLog is its mergeability. Two HLL instances of the same precision can be combined by taking the maximum value of each corresponding register, reflecting the union of elements from both data streams. This allows distributed or parallel systems to maintain separate sketches for subsets of data and aggregate them into a single unified estimate at any time. Another strength is the ability to handle any data type implementing Rust’s Hash trait, from simple strings to complex objects. Although HyperLogLog provides an approximate count rather than an exact value, its relative error typically remains low, and it can be tuned by adjusting the precision parameter p to trade off memory usage against accuracy.

The HyperLogLog structure is defined with a precision field p, a vector of 8-bit registers, and a constant alpha used in the cardinality estimation formula. During construction in new, the registers are initialized to zero, and alpha is set based on the number of registers. The add method takes an element, hashes it to 64 bits, and uses the top p bits to select a register. It shifts away these top bits and counts the leading zeros in the remainder to obtain a run-of-zeros value. If this value exceeds the register’s current value, the register is updated to reflect the newly observed maximum.

The estimation of unique elements happens in the estimate function, which applies the canonical HyperLogLog formula. The formula multiplies a bias-correcting term alpha by the square of the number of registers m, then divides by the sum of 2^(-register_value) across all registers. Because each register only holds the largest run of zeros encountered, this set of counters succinctly summarizes the distribution of hash values in the data. The from_iter constructor provides a convenient way to build a HyperLogLog by passing in an iterator of elements, while the merge method allows two sketches of the same precision to be combined, an operation that is particularly useful in large-scale or distributed applications. Finally, the main function demonstrates how to apply HyperLogLog to integer values, string data with duplicates, and how merging two separate sketches approximates the count of their combined sets.

#### Project Structure:

```plaintext
experiment_51_1/
├── Cargo.toml                  # Rust dependencies
└── src/
    ├── main.rs                 # Rust implementation
    └── output.txt              # Output file
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
only use std
```

#### Explanation of the Output

##### 1. Integer Count Estimation

```rust
Actual integer count: 10000
Estimated integer count: 9560.80
```

* A dataset containing 10,000 unique integers (0 to 9,999) is processed using HyperLogLog with precision p=10.
* The estimated count is 9560.80, which is close but slightly under the actual count due to the probabilistic nature of HyperLogLog.
* The estimation error is (10000 - 9560.80) / 10000 = 4.39%, which is expected as HyperLogLog is an approximation algorithm.

##### 2. String Count Estimation

```rust
Actual string unique count: 5
Estimated string count: 13.51
```

* A list of strings is processed: ["apple", "banana", "cherry", "banana", "date", "apple", "elderberry"].

* The actual number of unique strings is 5 ({"apple", "banana", "cherry", "date", "elderberry"}).

* The estimated count is 13.51, which is significantly higher than the actual value.

* The overestimation is likely due to the low precision (p=4), meaning there are only 16 registers (2^4 = 16), leading to more collisions and increased variance.

##### 3. Merging Two HyperLogLogs

```rust
Merging two HyperLogLogs each containing half of the integer range:
Merged estimate of unique integers: 9560.80
Actual unique count (0..10000): 10000
```

* The dataset of 10,000 integers is split into two halves (0..5000 and 5000..9999).

* Two separate HyperLogLogs are created for each half.

* After merging them, the estimated count remains 9560.80, the same as the previous estimate.

* This confirms that merging two HyperLogLogs does not overcount, as they retain only the maximum register values.

#### Conclusion

* HyperLogLog provides an efficient way to estimate cardinality (unique elements) using limited memory.

* The estimation is close but not exact—it depends on the precision (p value) and the nature of the data.

* Higher precision (p) reduces error but requires more memory.

* Smaller datasets (like strings with p=4) suffer from higher variance, leading to more overestimation or underestimation.

* Merging HyperLogLogs maintains accuracy, showing its usefulness for distributed data aggregation.

For better accuracy in real-world applications, choosing an appropriate p value is crucial, balancing memory usage and estimation precision.

