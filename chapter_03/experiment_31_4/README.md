## 3.1. Introduction to Data Structures and Algorithms

### experiment_31_4

This code demonstrates how to conditionally offload computations to GPU or FPGA in Rust by relying on feature flags in Cargo. In a single program, it provides three possible paths: one for GPU-accelerated processing, one for FPGA-accelerated processing, and a CPU-based fallback. Depending on which feature is enabled at compile time, the code will either initialize the relevant device and perform a placeholder transformation on an array of floats, or default to a simple CPU calculation if neither specialized device is available. This allows for a flexible design in which a single source tree can adapt to a variety of hardware environments without complex branching in the codebase.

In this code, the GPU and FPGA modules each simulate device initialization before applying a trivial computation to every element in the input data. The GPU path doubles each value, while the FPGA path adds ten. In practice, these stubs would be replaced by real kernels or bitstream invocations, with memory transfers to and from the device. If neither feature is present, the program multiplies the data by 1.1 on the CPU to provide a fallback. At runtime, Rust checks whether the feature flags “gpu” or “fpga” are active, conditionally compiling in the corresponding code blocks. This approach keeps the code simple yet extensible, allowing developers to maintain specialized paths for different acceleration options.

#### Files contents:
* experiment_31_2/Cargo.toml (Cargo.toml file for dependencies)
* src/
  * main.rs (rust script)
  * output.txt (output file)

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### [dependencies]

```toml
rand = "0.9.0"
rayon = "1.10.0"
```

#### [features]
```toml
gpu = []
fpga = []
```

##### Explanation of the Output:
This Rust program is designed to demonstrate how computations can be offloaded to different hardware accelerators (GPU, FPGA) or fallback to CPU processing, depending on the active feature flags at compile time. The program uses conditional compilation with cfg attributes to control which hardware-accelerated path gets executed.

##### Key Components:

###### 1. Feature Flags (gpu, fpga):

* The program supports three execution paths:
  * GPU acceleration (if the "gpu" feature is enabled)
  * FPGA acceleration (if the "fpga" feature is enabled)
  * CPU fallback (if neither "gpu" nor "fpga" is enabled)

###### 2. GPU Path (enabled via "gpu" feature):

* The gpu_impl module initializes the GPU environment and processes data using a placeholder GPU function (multiplying each value by 2.0).
* In a real-world scenario, this would involve actual GPU computations through libraries like CUDA or OpenCL.

###### 3. FPGA Path (enabled via "fpga" feature):

* The fpga_impl module simulates FPGA initialization and processing by adding 10.0 to each input value.
* In a real-world scenario, FPGA-specific calls would be made to offload computation onto FPGA hardware.

###### 4. CPU Path (when neither "gpu" nor "fpga" features are enabled):

* If neither GPU nor FPGA features are enabled, the program falls back to a simple CPU-based computation.
* The CPU fallback simulates a computation by multiplying each input value by 1.1. This computation is sequential unless you choose to parallelize it with something like Rayon (though it's not used here).

##### Program Flow:
* The program is executed with no GPU/FPGA features enabled, so the program falls back to the CPU path.
* The input data vector [1.0, 2.0, 3.0, 4.0, 5.0] is processed by the CPU path.
* Each value in the input data is multiplied by 1.1, resulting in:
  * 1.0 * 1.1 = 1.1
  * 2.0 * 1.1 = 2.2
  * 3.0 * 1.1 = 3.3
  * 4.0 * 1.1 = 4.4
  * 5.0 * 1.1 = 5.5
* The result is printed as:
```rust
CPU-only fallback result: [1.1, 2.2, 3.3000002, 4.4, 5.5]
```

##### Output Breakdown:
* "No GPU/FPGA features enabled. Using CPU fallback in parallel if desired."
  * This message is printed because neither the "gpu" nor "fpga" feature flags are enabled during compilation. It informs the user that the CPU fallback path is being used.
* CPU-only fallback result: [1.1, 2.2, 3.3000002, 4.4, 5.5]
  * This is the result after applying the CPU fallback function, where each value in the input data is multiplied by 1.1. The result shows slight floating-point precision variations (e.g., 3.3000002 instead of 3.3), which are common when working with floating-point arithmetic.

#### Conclusion:
* The program demonstrates a simple way to handle hardware acceleration using conditional compilation in Rust based on feature flags.
  * GPU and FPGA paths would be used if the respective features (gpu or fpga) were enabled during compilation. The program provides placeholder code for both paths.
  * CPU fallback is used when neither GPU nor FPGA features are active. This is useful for environments where no hardware acceleration is available or for benchmarking purposes.
* The CPU path performs a simple operation (multiplying by 1.1), but in real-world applications, more complex computations could be offloaded to the GPU or FPGA for better performance.
* Key takeaway: Depending on the compile-time feature flags, the program can either offload computation to specialized hardware (GPU/FPGA) or fall back to CPU processing. This approach allows developers to write efficient, portable code that adapts to different hardware configurations.



