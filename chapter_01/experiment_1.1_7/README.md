## 2.1. Introduction to Rust Programming Language

### experiment_21_7

Here is a brief example illustrating how Rust can handle an imperative, pointer-based approach to building a simple suffix array. This style is reminiscent of C++ in that it manipulates raw pointers for performance-critical tasks. Rust’s compiler still enforces strict rules about the pointer’s lifetime, preventing you from using it once the underlying data goes out of scope.

In this snippet, build_suffix_array manually handles pointer arithmetic to compare suffixes, reminiscent of a lower-level C++ style. The reference genome data is stored in a Vec, and seq.as_ptr() provides a raw pointer to the first byte. Because Rust attaches the pointer’s lifetime to the underlying Vec, any misuse—such as accessing memory after seq goes out of scope—would result in a compile-time error. The compare_suffixes function uses unsafe code, but the scope of that unsafety is minimized. Rust checks that the pointer remains valid as long as the Vec has not been dropped, preserving memory safety guarantees even though we resort to pointer-based comparisons for performance or fine-grained control.

#### Project Structure:

```plaintext
experiment_21_7/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── reads.fasta                # FASTA file containing sequence reads
    └── output.txt                 # Output file
```

#### How to run:

```rust
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

only use "std"

#### Explanation of the Output

##### Understanding the Code
This Rust program constructs a suffix array for a given reference genome sequence (ATGCGT). A suffix array is a sorted list of indices representing all possible suffixes of a string, sorted in lexicographical (dictionary) order.

##### 1. Reference Genome

```rust
let reference = b"ATGCGT".to_vec();
```

* This is a DNA-like sequence stored as a Vec<u8> (b"ATGCGT" represents bytes for "ATGCGT").
* We generate suffixes from this sequence and sort them lexicographically.

##### 2. Building the Suffix Array

```rust
let suffix_array = build_suffix_array(&reference);
```

* This function constructs the suffix array by sorting suffixes of the string based on their lexicographical order.
* It utilizes raw pointers for memory-efficient substring comparison.

##### 3. Sorting Suffixes Using compare_suffixes

```rust
unsafe fn compare_suffixes(ptr: *const u8, len: usize, i: usize, j: usize) -> Ordering
``

* This function compares two suffixes by directly accessing memory using raw pointers.
* The cmp() function is used to determine their lexicographical order.

##### Breaking Down the Output

For the input string "ATGCGT", here are all suffixes with their starting positions:

```sh
Index	Suffix
0	    ATGCGT
1	    TGCGT
2	    GCGT
3	    CGT
4	    GT
5	    T
```

When sorted lexicographically:

1. ATGCGT (Index 0)
2. CGT (Index 3)
3. GCGT (Index 2)
4. GT (Index 4)
5. T (Index 5)
6. TGCGT (Index 1)
7. 
Thus, the suffix array output:

```rust
Suffix Array for "ATGCGT": [0, 3, 2, 4, 5, 1]
```

This means:
* The 0th character (ATGCGT) is lexicographically smallest.
* The 3rd character (CGT) is the next smallest.
* The 2nd character (GCGT) follows, and so on.
  
#### Conclusion
* The suffix array correctly sorts all suffixes lexicographically and outputs their starting indices.
* This approach, although naive, demonstrates the fundamental suffix array construction.
* Raw pointers (unsafe code) allow efficient memory access for direct substring comparison.
* Suffix arrays are widely used in bioinformatics (genome indexing, sequence alignment) and string matching algorithms (like Burrows-Wheeler Transform, FM-index).

This implementation can be optimized further using advanced suffix array construction algorithms like Kärkkäinen-Sanders, SA-IS, or suffix trees.










