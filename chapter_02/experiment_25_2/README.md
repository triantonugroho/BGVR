## 2.5. Acquiring and Cleaning Data

### experiment_25_2

In this pipeline, all software is assumed to be installed on your machine or HPC cluster. If you have a module system on HPC, you might load specific Rust and Samtools versions by name. On a local PC, you might rely on a version manager or direct binary installation. Either way, the principle remains the same: you explicitly declare which versions of your tools you need so that colleagues or future runs of the pipeline can replicate the environment as closely as possible. By tracking data version labels, numeric parameters, and environment variables, you create a clear, traceable record of how each result was produced.

Below is a single Nextflow script that illustrates these concepts. You can store this script in a file named main.nf.

In this script, both the Rust and Samtools versions are pinned as parameters, but the actual mechanics of ensuring those exact versions are used depends on how your local or HPC system is managed. If you are on an HPC cluster with a module system, you might insert commands like module load rust/1.66.0 or module load samtools/1.12 within the script blocks, or you might define them in a nextflow.config file as directives for each process. If you are on a local workstation, you would install Rust 1.66.0 (for example, via rustup) and Samtools 1.12 (via a package manager or manual build) to match these pinned versions. In a fully production setting, you would document precisely how to install or load these versions so that any collaborator can do so in the future.

The script also illustrates data versioning by labeling the input dataset as dataset_v1.0 and passing that label as params.input_data. In a real pipeline, you could store data files under versioned directories or attach metadata that matches this label, ensuring that your logs or outputs clearly indicate which dataset was used. Additionally, the example sets a numeric parameter (params.numeric_param) to 21, which might represent a hyperparameter for a Rust-based algorithm or Samtools filter threshold. Environment variables are also demonstrated by printing and exporting params.env_variable so that the Rust build step or subsequent tasks can adjust their behavior accordingly.

Each step in Nextflow is a process that consumes and produces data via channels. In the BUILD_RUST_APP process, the pipeline compiles a Rust program, storing the resulting binary in myapp and publishing it to a publish_output directory. When the pipeline transitions to the RUN_TOOL process, the newly compiled binary is provided as input, and the script verifies the Samtools version before running the binary. Outputs produced during this step, in this case results.txt, are also published to publish_output. The final workflow block collects the resulting files but can be extended with additional processes if needed.

#### Files contents:
* main.nf (nextflow script)
* src/dataset_v1.0 (dataset folder)
* src/publish_ouput/
  * myapp
  * results.txt
* nextflow.log.7 (nextflow log file)
* Cargo.toml (Cargo.toml file)
* output.txt (output file)

#### How to run:

run in WSL:

```wsl
nextflow run main.nf | tee output.txt
```

(run main.nf and save the output in output.txt)
  
#### [dependencies]

no dependencies

#### Explanation of the Output

This Nextflow workflow (main.nf) automates the process of:
1. Compiling a Rust application using a specific Rust version.
2. Running the compiled Rust application alongside Samtools while using a numeric parameter.

##### Step-by-Step Execution

Each step in the workflow runs as an independent Nextflow process, and their execution is tracked.

###### 1. Workflow Initialization

```nextflow
N E X T F L O W   ~  version 24.10.4
Launching `main.nf` [prickly_ramanujan] DSL2 - revision: 330781e9fa
```

* Nextflow version: 24.10.4
* Workflow identifier: prickly_ramanujan (Nextflow assigns a random name to each execution).
* DSL2 enabled: The workflow is structured using modular pipeline design.
  
###### 2. Compiling the Rust Application

```nextflow
executor >  local (1)
[c8/908b38] BUILD_RUST_APP [  0%] 0 of 1
```

* The BUILD_RUST_APP process starts.
* Initially, 0 of 1 Rust applications have been built.

```nextflow
executor >  local (1)
[c8/908b38] BUILD_RUST_APP [100%] 1 of 1 ✔
```

* BUILD_RUST_APP completed successfully (100%).
* The Rust compiler (rustc) created an executable named myapp.
* This file is now ready for execution in the next process.

###### 3. Running the Compiled Rust Application

```nextflow
executor >  local (2)
[c8/908b38] BUILD_RUST_APP [100%] 1 of 1 ✔
[62/7bc76e] RUN_TOOL       [  0%] 0 of 1
```

* RUN_TOOL process starts.
* The Rust application myapp and the numeric parameter are passed as input.

```nextflow
executor >  local (2)
[c8/908b38] BUILD_RUST_APP [100%] 1 of 1 ✔
[62/7bc76e] RUN_TOOL       [100%] 1 of 1 ✔
```

* The Rust application executed successfully (100%).
* The Rust program printed "Hello from pinned Rust!" as expected.
* The numeric parameter (21) was also used in the output:

```nextflow
Results for param = 21
```

* The results were saved in results.txt.

###### 4. Workflow Completion Summary

```nextflow
executor >  local (2)
[c8/908b38] BUILD_RUST_APP [100%] 1 of 1 ✔
[62/7bc76e] RUN_TOOL       [100%] 1 of 1 ✔
```

* All processes completed successfully (100%).

#### Conclusion

* The workflow successfully completed all steps:
  * Compiled a Rust application using Rust version 1.85.0.
  * Executed the Rust application while using a numeric parameter (21).
  * Samtools version 1.19.2 was verified, ensuring reproducibility.
  * Final results were saved in results.txt.
* Key Takeaways:
  * Reproducible environment: Rust and Samtools versions are explicitly set.
  * Modular design: Each step is independent and can be modified easily.
  * Workflow correctness: The expected outputs were generated, confirming a successful run. 









