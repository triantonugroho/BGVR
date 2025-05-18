## 2.1. Introduction to Rust Programming Language

### experiment_21_9

Here is a simple example demonstrating how a Rust library can use tch (the Rust bindings for PyTorch) and then expose this functionality to Python via PyO3, allowing Rust to handle computationally heavy workloads while Python remains the interface for scripting or orchestrating workflows. First, create a new Cargo project and edit the Cargo.toml file so that the library name matches pyo3_tch_example, specifying cdylib for the crate type. Include both PyO3 and tch in the dependencies, enabling the extension-module feature for PyO3 so that a Python extension can be compiled. In src/lib.rs, place the following code.

This code defines a Python module called pyo3_tch_example, which includes a single function run_linear_forward_py. The function creates a linear layer in Rust, transforms the input vector into a one-by-two tensor on the CPU, applies a forward pass, and returns the result as a Vec. You can build this project with either cargo build --release or maturin build --release. Once built, you will obtain a compiled shared library (for example, pyo3_tch_example.so on Linux) that can be imported directly into Python:

```python
import pyo3_tch_example as tch_ex

# Pass exactly two floats for the forward pass.
result = tch_ex.run_linear_forward_py([1.0, -0.5])
print("Rust + tch output:", result)
```

This setup demonstrates how PyO3 can bridge Rust and Python, letting you leverage Rust’s concurrency and memory safety guarantees for computationally heavy tasks while seamlessly integrating with Python-based machine learning or data analysis pipelines.

#### Project Structure:

```plaintext
experiment_21_9/
├── lib.rs                         # Rust script to make Python library
├── Testing pyo3_tch_example.ipynb # Python notebook for testing
├── Cargo.toml                     # Rust project configuration and dependencies
├── target/debug/
│   ├── pyo3_tch_example.d
│   ├── pyo3_tch_example.dll
│   ├── pyo3_tch_example.dll.exp
│   ├── pyo3_tch_example.dll.lib
│   └── pyo3_tch_example.pdb
└── target/wheels/
    └── experiment_21_9-0.1.0-cp311-cp311-win_amd64.whl

~anaconda3/Lib/site-packages/pyo3_tch_example/
├── __init__                       # Python Source File
├── pyo3_tch_example.cp311-win_amd64.pyd
└── __pycache__/
    └── __init__.cpython-311      # Compiled Python File
```

#### How to run:

```powershell
cargo run 
```

(run lib.rs and make pyo3_tch_example python library)
  
#### [dependencies]

```toml
pyo3 = { version = "0.23.5", features = ["extension-module"] }
tch = "0.19.0"
```
#### [lib]

```toml
name = "pyo3_tch_example"
crate-type = ["cdylib"]
```
#### Explanation of the Output

##### Understanding the Code
This Rust program, using PyO3 and tch (Torch for Rust), creates a Python library named pyo3_tch_example. The library provides a function (run_linear_forward_py) that applies a simple linear transformation (y = Wx + b) using PyTorch in Rust.

##### Key Code Components

##### 1. Function: run_linear_forward_py

* Accepts a vector of two float values (Vec<f32>) as input.
* Checks if the input contains exactly two elements (error otherwise).
* Uses PyTorch (tch) to create a linear layer with 2 input neurons and 1 output neuron.
* Converts the input into a 1×2 PyTorch tensor and runs a forward pass.
* Returns the output tensor as a Rust vector (Vec<f32>).

##### 2. Creating a Linear Layer in PyTorch (tch)

```rust
let linear = nn::linear(&root, 2, 1, Default::default());
```

* Defines a fully connected layer with 2 input neurons and 1 output neuron.
* PyTorch automatically initializes weights (W) and bias (b).

##### 5. Processing Input Data

```rust
let input_tensor = Tensor::from_slice(&input_data)
     .view([1, 2]) // Reshape into a 1x2 tensor
     .to_device(device); // Use CPU
```

* Converts the input (Vec<f32>) into a 1×2 PyTorch tensor.
* Uses .view([1, 2]) to reshape the tensor.

##### 4. Running the Forward Pass

```rust
let output = linear.forward(&input_tensor);
```

* Computes y = Wx + b using the initialized weights and bias.

##### 5. Converting Output Back to Rust Vector

```rust
let out_vec: Vec<f32> = output.try_into().map_err(|e| {
    PyErr::new::<PyValueError, _>(format!("Tensor conversion error: {}", e))
})?;
```

* Extracts the output tensor and converts it into a Rust vector (Vec<f32>).
  
##### Expected Output
The function is registered as a Python module, allowing Python users to call:

```python
import pyo3_tch_example
output = pyo3_tch_example.run_linear_forward_py([1.0, 2.0])
print(output)
```

* The output is a list with a single float value (e.g., [0.345]), which depends on PyTorch's random initialization of weights and bias.
 
#### Conclusion
* The pyo3_tch_example library allows seamless integration of Rust and PyTorch into Python.
* The function performs a simple linear transformation (like a single neuron with 2 inputs).
* This approach enables fast, efficient deep learning computations using Rust while maintaining Python accessibility.
* The function can be extended to support more complex neural network operations.






