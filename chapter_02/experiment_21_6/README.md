## 2.1. Introduction to Rust Programming Language

### experiment_21_6

Below is a simple example that demonstrates Rust’s zero-cost abstractions in a genomic setting. This example reads a FASTA file and uses Rust’s iterator adapters (such as map, filter, and for_each) to process genomic sequences without incurring additional runtime overhead. Despite the high-level functional style, the compiler optimizes these operations down to efficient machine code, making them comparable to hand-written loops in languages like C or C++.

This code reads a FASTA file using the bio::io::fasta crate, which emits a stream of records (each containing a sequence). For each record, the sequence bytes are converted into a String via String::from_utf8_lossy, then any sequence under 50 nucleotides is discarded through a filter call. The remaining sequences move to a second map step that calculates their GC content by counting the characters ‘G’ or ‘C.’ Finally, the .sum() operation aggregates these individual GC counts into one total. Rust compiles this chain of iterator adapters into efficient loops without constructing intermediary collections, a technique referred to as “zero-cost abstraction.” Consequently, although the code is written in a high-level, functional style, it executes at speeds comparable to hand-tuned loops in lower-level languages.

#### Files contents:
* lib.rs (rust script to make library python)
* Testing pyo3_biopython_example.ipynb (python notebook)
* Cargo.toml (Cargo.toml file)
* target/debug/
  * pyo3_biopython_example.d
  * pyo3_biopython_example.dll
  * pyo3_biopython_example.dll.exp
  * pyo3_biopython_example.dll.lib
  * pyo3_biopython_example.pdb
* target/wheels/
  * experiment_21_6-0.1.0-cp311-cp311-win_amd64.whl
* ~anaconda3/Lib/site-packages/pyo3_biopython_example/
  * __init__ (Python Source File)
  * pyo3_biopython_example.cp311-win_amd64.pyd
* ~anaconda3/Lib/site-packages/pyo3_biopython_example/__pycache__/
  * __init__.cpython-311 (Compiled Python File)

#### How to run:

```powershell
cargo run 
```

(run lib.rs and make pyo3_biopython_example python library)
  
#### [dependencies]

```toml
pyo3 = { version = "0.23.5", features = ["extension-module"] }
```
#### [lib]

```toml
name = "pyo3_biopython_example"
crate-type = ["cdylib"]
```
#### Explanation of the Output
The output "pyo3_biopython_example (python library)" indicates that the Rust code has been successfully compiled into a Python library using PyO3. This means that Python users can now import and use the Rust functions as if they were native Python functions.

##### Breakdown of the Rust Code and Its Purpose

This Rust library defines two bioinformatics-related functions:

```sh
reverse_complement(seq: &str) -> PyResult<String>
```

* Takes a DNA sequence (seq) as input.
* Computes the reverse complement (i.e., replaces each nucleotide with its complementary base and reverses the sequence).
* Example:

```sh
>>> reverse_complement("ATGC")
"GCAT"
```

```sh
gc_content(seq: &str) -> PyResult<f64>
```

* Calculates the GC content percentage of a given DNA sequence.
* GC content is important in genetics because it influences DNA stability.
* Example:

```sh
>>> gc_content("ATGC")
50.0
```

* If an empty sequence is provided, it raises a ValueError.

```sh
pyo3_biopython_example(py: Python, m: &Bound<', PyModule>) -> PyResult<()>
```

* This function creates a Python module named pyo3_biopython_example.
* It exposes the two Rust functions (reverse_complement and gc_content) to Python.
* When compiled, it can be imported in Python as:

```python
import pyo3_biopython_example
print(pyo3_biopython_example.reverse_complement("ATGC"))  # Output: "GCAT"
print(pyo3_biopython_example.gc_content("ATGC"))  # Output: 50.0
```

#### Conclusion
The Rust code successfully integrates with Python via PyO3, making it usable as a Python module.
This allows users to execute high-performance Rust functions in Python for bioinformatics tasks like reverse complement computation and GC content analysis.
The combination of Rust’s speed and Python’s ease of use makes this approach ideal for large-scale biological sequence analysis.









