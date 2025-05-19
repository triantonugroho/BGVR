## 1.1. Introduction to Rust Programming Language

### experiment_1_1_6

When building a Rust project for Python interoperability, you can specify crate-type = ["cdylib"] to ensure your library is compiled as a dynamic library suitable for Python imports. Enabling PyO3 with the "extension-module" feature allows you to generate a Python C-extension, making your library compatible with the Python ecosystem. Additionally, you may include maturin as a build dependency to simplify packaging your Rust library into Python wheels, streamlining distribution and installation.

In this code, you configure a Rust crate as a Python extension by setting crate-type = ["cdylib"] and enabling PyO3’s "extension-module" feature. In reverse_complement, the snippet imports BioPython’s Bio.Seq module, dynamically creates a Seq object from the Rust string, and calls the reverse_complement method. Any Python exceptions become PyErrs, gracefully bridging the Python–Rust boundary. The second function, gc_content, shows how you can call additional modules (e.g., Bio.SeqUtils.GC) in exactly the same manner. This design is robust (due to explicit error handling), comprehensive (extensible by simply defining new Rust functions and exposing them to Python), and scalable (capable of multi-threading in Rust for CPU-bound tasks and seamlessly returning control to Python for broader pipeline logic). Once compiled using Maturin or Cargo, you can import the .so (or .pyd) file in Python, meaning you get a fully integrated environment combining Python’s extensive libraries with Rust’s performance.

#### Project Structure:

```plaintext
experiment_1_1_6/
├── lib.rs                         # Rust script to make Python library
├── Testing pyo3_biopython_example.ipynb # Python notebook for testing
├── Cargo.toml                     # Rust project configuration and dependencies
├── target/debug/
│   ├── pyo3_biopython_example.d
│   ├── pyo3_biopython_example.dll
│   ├── pyo3_biopython_example.dll.exp
│   ├── pyo3_biopython_example.dll.lib
│   └── pyo3_biopython_example.pdb
└── target/wheels/
    └── experiment_1_1_6-0.1.0-cp311-cp311-win_amd64.whl

~anaconda3/Lib/site-packages/pyo3_biopython_example/
├── __init__                       # Python Source File
├── pyo3_biopython_example.cp311-win_amd64.pyd
└── __pycache__/
    └── __init__.cpython-311      # Compiled Python File
```

#### Cargo.toml

```toml
[package]
name = "experiment_1_1_6"
version = "0.1.0"
edition = "2021"

[dependencies]
pyo3 = { version = "0.23.5", features = ["extension-module"] }

[lib]
name = "pyo3_biopython_example"
crate-type = ["cdylib"]
```

#### How to run:

```powershell
cargo run 
```

(run lib.rs and make pyo3_biopython_example python library)
  
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









