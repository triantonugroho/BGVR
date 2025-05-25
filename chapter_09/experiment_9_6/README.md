## 9.6 Advanced Topics: Single-Cell RNA-seq and Beyond

### experiment_9_6

Below is a Rust snippet illustrating how one might load a sparse matrix of single-cell gene counts, apply a simple dimensionality reduction (via truncated PCA), and then export the resulting cell coordinates. This example uses “sprs” for sparse matrices, “ndarray-linalg” for linear algebra routines, and “serde” for data input/output.

This code first reads a sparse count file, constructs a “sprs” sparse matrix, and then converts it to a dense matrix for demonstration purposes. Because single-cell datasets can be extremely large, true production workflows typically employ more advanced, memory-efficient methods for SVD, such as randomized or iterative approaches. The sprs crate provides a variety of sparse matrix utilities, and ndarray-linalg facilitates linear algebra routines on ndarray structures. For industrial-scale usage, it is prudent to incorporate concurrency (via rayon) and robust error handling.

Below is a Nextflow script showcasing how to incorporate this Rust step into a broader scRNA-seq pipeline that might include tasks such as demultiplexing, alignment, and the generation of a sparse count matrix.

Such a pipeline highlights the modular design that is typical in Nextflow, where each step performs a distinct function, and the final output can be used for downstream visualization or clustering. Larger HPC clusters or cloud-based solutions commonly harness Nextflow’s inherent scalability, automatically distributing tasks to run in parallel or at scale.

Many AI engineers and bioinformaticians now integrate single-cell pipelines within Nextflow to manage immense scRNA-seq datasets, particularly when investigating rare cell populations relevant to complex diseases. One global consortium found that porting core data transformations to Rust significantly reduced execution times for large single-cell projects, enabling them to analyze millions of cells with minimal downtime. This streamlined approach accelerated iterative research cycles, ultimately helping identify novel cell subtypes implicated in inflammatory pathologies. By pairing the computational efficiency of Rust with Nextflow’s workflow management, teams can tackle the unique challenges of single-cell data at scale.

#### Project Structure:

```plaintext
experiment_9_6/
├── Cargo.toml                               # Rust package configuration and dependencies
├── main.nf                                  # Nextflow pipeline script
├── test_analysis.txt                        # Test analysis output
├── test-normalized.tsv                      # Test normalization output
├── test_summary                             # Test summary output
├── README.md                                # Project documentation
│
├── data/                                    # Generated dataset folder
|   ├── annotation/                          # Annotation dataset folder
|   │   └── annotation.gtf                   # Annotation data
|   ├── genome_index/                        # Genome index dataset folder 
|   │   ├── chrNameLength.txt                # Character name length data
|   │   └── genomeParameters.txt             # Genome parameters data
|   ├── reads/                               # Reads dataset folder 
|   |    └── sample1_1.fastq.gz              # sample 1_1 fastq file
|   |    └── sample1_2.fastq.gz              # sample 1_2 fastq file
|   |    └── sample2_1.fastq.gz              # sample 2_1 fastq file
|   |   ├── sample2_1.fastq.gz               # sample 2_2 fastq file
|   |   ├── sample3_1.fastq.gz               # sample 3_1 fastq file
|   |   ├── sample3_1.fastq.gz               # sample 3_2 fastq file
|   |   ├── sample4_1.fastq.gz               # sample 4_1 fastq file
|   |   └── sample4_1.fastq.gz               # sample 4_2 fastq file
|   ├── expression_data.tsv                  # Expression data
|   ├── gene_info.tsv                        # Gene information data
|   └── gene_info.tsv                        # Gene information data
├── results/                                 # Results/output folder
|   ├── expression/                          # Expression result folder
|   │   ├── sample1_expression.tsv           # Sample 1 expression result
|   │   ├── sample2_expression.tsv           # Sample 2 expression result
|   │   ├── sample3_expression.tsv           # Sample 3 expression result
|   │   └── sample4_expression.tsv           # Sample 4 expression result
|   ├── final/                               # Final result folder 
|   │   ├── analysis_summary.txt             # Analysis summary result
|   │   └── combined_matrix.tsv              # Combined matrix result
|   ├── qc/                                  # Reads dataset folder 
|   |   ├── sample1_qc.txt                   # Sample 1 qc result
|   |   ├── sample2_qc.txt                   # Sample 2 qc result
|   |   ├── sample3_qc.txt                   # Sample 3 qc result
|   |   └── sample4_qc.txt                   # Sample 4 qc result
|   └── reports/                             # Reports folder 
|       ├── execution_report.html            # HTML execution report
|       ├── execution_timeline.html          # HTML execution timeline
|       └── execution_trace.txt              # Execution trace 
├── src/                                     # Rust source code
│   └── main.rs                              # Main Rust expression tool implementation
│
├── target/                                  # Rust build artifacts
│   └── release/
│       └── rust_expression_tool             # Compiled Rust executable binary
│
└── work/                                    # Nextflow working directory (temporary files)
    ├── 0e/
    │   └── 7e2ae70a91328dfb08b6055a18734b/
    │       └── sample1_expression           # Sample 1 expression result
    ├── 17/
    │   └── 87b14f09d1f0ef65cd63bbb72e30d6/
    │       └── sample4_qc.txt               # Sample 4 qc result
    ├── 1e/
    │   └── 00fee1ccf630ae1132ff0fa89e271d/
    │       └── sample1_qc.txt               # Sample 1 qc result
    ├── 24/
    │   └── 5f385617f7771f0a036e582039c377/
    │       └── sample3_expression.tsv       # Sample 3 expression result
    ├── 2f/
    │   └── 069782c3fb1fb1d8ce94e306d9c51a/
    │       ├── analysis_summary.txt         # Analysis summary result
    │       ├── combined_matrix.tsv          # Combined matrix result
    │       └── sample_metadata.tsv          # Sample metadata result
    │   └── 00fee1ccf630ae1132ff0fa89e271d/
    │       └── sample4_expression.tsv       # Sample 4 expression result
    ├── 32/
    │   └── cf7f104a40c75f28daeb804176a4b6/
    │       └── sample1_expression.tsv       # Sample 1 expression result
    ├── 49/
    │   └── ba2ca3654969d439a2d7c2cd59eac8/
    │       └── sample2_expression.tsv       # Sample 2 expression result
    ├── 4f/
    │   └── f123bdd7ae27d534a2ba36d86cc01b/
    │       └── sample1_qc.txt               # Sample 1 qc result
    ├── 6d/
    │   └── 7790e3f714c345dc2c73df352efa34/
    │       └── sample3_qc.txt               # Sample 3 qc result
    ├── 79/
    │   └── 5fa81ff2c0fdc2db9c83231db5cf76/
    │       └── sample3_expression.tsv       # Sample 3 expression result
    ├── 88/
    │   └── 2b569525cb9550201e448b34bdab5f/
    │       └── sample3_qc.txt               # Sample 3 qc result
    ├── 9f/
    │   └── 09c10b61ab5dea1edfc43a70dff06b/
    │       ├── analysis_summary.txt         # Analysis summary result
    │       ├── combined_matrix.tsv          # Combined matrix result
    │       └── sample_metadata.tsv          # Sample metadata result    
    ├── bc/
    │   └── 7689d4ecab5b158eb7c73fcfe55f0f/
    │       └── sample2_qc.txt               # Sample 2 qc result
    ├── c5/
    │   └── 45063c17268c9a506d05248d560361/
    │       └── sample2_qc.txt               # Sample 2 qc result
    ├── e5/
    │   └── 698ae5eb2ff9a10e60e73ebd4578f6/
    │       └── sample4_expression.tsv       # Sample 4 expression result
    ├── fa/
    │   └── 92f7de84fa4b40997c34ea17fbdfe0/
    │       └── sample4_qc.txt               # Sample 4 qc result
    ├── fb/
        └── cc864e3f9b36df849f9799028f87f5/
            └── sample2_expression.tsv       # Sample 2 expression result
```

#### Cargo.toml

```toml
[package]
name = "differential-expression-analyzer"
version = "1.0.0"
edition = "2021"
authors = ["Bioinformatics Pipeline <pipeline@example.com>"]
description = "A robust differential expression analysis tool for RNA-seq data"
license = "MIT"
repository = "https://github.com/username/differential-expression-analyzer"
keywords = ["bioinformatics", "rnaseq", "differential-expression", "genomics", "statistics"]
categories = ["science", "command-line-utilities"]

[[bin]]
name = "diff-expr-analyzer"
path = "src/main.rs"

[dependencies]
clap = "4.4"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
ndarray = "0.15"
statrs = "0.16"
log = "0.4"
env_logger = "0.10"
csv = "1.3"
anyhow = "1.0"
thiserror = "1.0"

[dev-dependencies]
tempfile = "3.8"
assert_cmd = "2.0"
predicates = "3.0"
approx = "0.5"

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
panic = "abort"

[profile.dev]
opt-level = 0
debug = true
```

#### How to run:

##### Generate sample data in wsl:

```wsl

# Update pip first

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ python3 -m pip install --upgrade pip
Requirement already satisfied: pip in /home/trian/miniconda3/lib/python3.12/site-packages (25.0)
Collecting pip
  Downloading pip-25.1.1-py3-none-any.whl.metadata (3.6 kB)
Downloading pip-25.1.1-py3-none-any.whl (1.8 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1.8/1.8 MB 382.0 kB/s eta 0:00:00
Installing collected packages: pip
  Attempting uninstall: pip
    Found existing installation: pip 25.0
    Uninstalling pip-25.0:
      Successfully uninstalled pip-25.0
Successfully installed pip-25.1.1

# Install scientific computing packages

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ python3 -m pip install --user numpy pandas scipy matplotlib seaborn scikit-learn
Requirement already satisfied: numpy in /home/trian/miniconda3/lib/python3.12/site-packages (2.2.5)
Requirement already satisfied: pandas in /home/trian/.local/lib/python3.12/site-packages (2.2.3)
Collecting scipy
  Downloading scipy-1.15.3-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (61 kB)
Requirement already satisfied: matplotlib in /home/trian/.local/lib/python3.12/site-packages (3.10.3)
Requirement already satisfied: seaborn in /home/trian/.local/lib/python3.12/site-packages (0.13.2)
Collecting scikit-learn
  Downloading scikit_learn-1.6.1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (18 kB)
Requirement already satisfied: python-dateutil>=2.8.2 in /home/trian/.local/lib/python3.12/site-packages (from pandas) (2.9.0.post0)
Requirement already satisfied: pytz>=2020.1 in /home/trian/.local/lib/python3.12/site-packages (from pandas) (2025.2)
Requirement already satisfied: tzdata>=2022.7 in /home/trian/.local/lib/python3.12/site-packages (from pandas) (2025.2)
Requirement already satisfied: contourpy>=1.0.1 in /home/trian/.local/lib/python3.12/site-packages (from matplotlib) (1.3.2)
Requirement already satisfied: cycler>=0.10 in /home/trian/.local/lib/python3.12/site-packages (from matplotlib) (0.12.1)
Requirement already satisfied: fonttools>=4.22.0 in /home/trian/.local/lib/python3.12/site-packages (from matplotlib) (4.58.0)
Requirement already satisfied: kiwisolver>=1.3.1 in /home/trian/.local/lib/python3.12/site-packages (from matplotlib) (1.4.8)
Requirement already satisfied: packaging>=20.0 in /home/trian/miniconda3/lib/python3.12/site-packages (from matplotlib) (24.2)
Requirement already satisfied: pillow>=8 in /home/trian/.local/lib/python3.12/site-packages (from matplotlib) (11.2.1)
Requirement already satisfied: pyparsing>=2.3.1 in /home/trian/.local/lib/python3.12/site-packages (from matplotlib) (3.2.3)
Collecting joblib>=1.2.0 (from scikit-learn)
  Downloading joblib-1.5.1-py3-none-any.whl.metadata (5.6 kB)
Collecting threadpoolctl>=3.1.0 (from scikit-learn)
  Downloading threadpoolctl-3.6.0-py3-none-any.whl.metadata (13 kB)
Requirement already satisfied: six>=1.5 in /home/trian/.local/lib/python3.12/site-packages (from python-dateutil>=2.8.2->pandas) (1.17.0)
Downloading scipy-1.15.3-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (37.3 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 37.3/37.3 MB 335.9 kB/s eta 0:00:00
Downloading scikit_learn-1.6.1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (13.1 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 13.1/13.1 MB 364.3 kB/s eta 0:00:00
Downloading joblib-1.5.1-py3-none-any.whl (307 kB)
Downloading threadpoolctl-3.6.0-py3-none-any.whl (18 kB)
Installing collected packages: threadpoolctl, scipy, joblib, scikit-learn
Successfully installed joblib-1.5.1 scikit-learn-1.6.1 scipy-1.15.3 threadpoolctl-3.6.0

# Install additional useful packages

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ python3 -m pip install --user jupyter plotly
Collecting jupyter
  Downloading jupyter-1.1.1-py2.py3-none-any.whl.metadata (2.0 kB)
Collecting plotly
  Downloading plotly-6.1.1-py3-none-any.whl.metadata (6.9 kB)
Collecting notebook (from jupyter)
  Downloading notebook-7.4.2-py3-none-any.whl.metadata (10 kB)
Collecting jupyter-console (from jupyter)
  Downloading jupyter_console-6.6.3-py3-none-any.whl.metadata (5.8 kB)
Collecting nbconvert (from jupyter)
  Downloading nbconvert-7.16.6-py3-none-any.whl.metadata (8.5 kB)
Collecting ipykernel (from jupyter)
  Downloading ipykernel-6.29.5-py3-none-any.whl.metadata (6.3 kB)
Collecting ipywidgets (from jupyter)
  Downloading ipywidgets-8.1.7-py3-none-any.whl.metadata (2.4 kB)
Collecting jupyterlab (from jupyter)
  Downloading jupyterlab-4.4.2-py3-none-any.whl.metadata (16 kB)
Collecting narwhals>=1.15.1 (from plotly)
  Downloading narwhals-1.40.0-py3-none-any.whl.metadata (11 kB)
Requirement already satisfied: packaging in /home/trian/miniconda3/lib/python3.12/site-packages (from plotly) (24.2)
Collecting comm>=0.1.1 (from ipykernel->jupyter)
  Downloading comm-0.2.2-py3-none-any.whl.metadata (3.7 kB)
Collecting debugpy>=1.6.5 (from ipykernel->jupyter)
  Downloading debugpy-1.8.14-cp312-cp312-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (1.3 kB)
Collecting ipython>=7.23.1 (from ipykernel->jupyter)
  Downloading ipython-9.2.0-py3-none-any.whl.metadata (4.4 kB)
Collecting jupyter-client>=6.1.12 (from ipykernel->jupyter)
  Downloading jupyter_client-8.6.3-py3-none-any.whl.metadata (8.3 kB)
Collecting jupyter-core!=5.0.*,>=4.12 (from ipykernel->jupyter)
  Downloading jupyter_core-5.7.2-py3-none-any.whl.metadata (3.4 kB)
Collecting matplotlib-inline>=0.1 (from ipykernel->jupyter)
  Downloading matplotlib_inline-0.1.7-py3-none-any.whl.metadata (3.9 kB)
Collecting nest-asyncio (from ipykernel->jupyter)
  Downloading nest_asyncio-1.6.0-py3-none-any.whl.metadata (2.8 kB)
Collecting psutil (from ipykernel->jupyter)
  Downloading psutil-7.0.0-cp36-abi3-manylinux_2_12_x86_64.manylinux2010_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (22 kB)
Collecting pyzmq>=24 (from ipykernel->jupyter)
  Downloading pyzmq-26.4.0-cp312-cp312-manylinux_2_28_x86_64.whl.metadata (6.0 kB)
Collecting tornado>=6.1 (from ipykernel->jupyter)
  Downloading tornado-6.5.1-cp39-abi3-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (2.8 kB)
Collecting traitlets>=5.4.0 (from ipykernel->jupyter)
  Downloading traitlets-5.14.3-py3-none-any.whl.metadata (10 kB)
Collecting decorator (from ipython>=7.23.1->ipykernel->jupyter)
  Downloading decorator-5.2.1-py3-none-any.whl.metadata (3.9 kB)
Collecting ipython-pygments-lexers (from ipython>=7.23.1->ipykernel->jupyter)
  Downloading ipython_pygments_lexers-1.1.1-py3-none-any.whl.metadata (1.1 kB)
Collecting jedi>=0.16 (from ipython>=7.23.1->ipykernel->jupyter)
  Downloading jedi-0.19.2-py2.py3-none-any.whl.metadata (22 kB)
Collecting pexpect>4.3 (from ipython>=7.23.1->ipykernel->jupyter)
  Downloading pexpect-4.9.0-py2.py3-none-any.whl.metadata (2.5 kB)
Collecting prompt_toolkit<3.1.0,>=3.0.41 (from ipython>=7.23.1->ipykernel->jupyter)
  Downloading prompt_toolkit-3.0.51-py3-none-any.whl.metadata (6.4 kB)
Requirement already satisfied: pygments>=2.4.0 in /home/trian/miniconda3/lib/python3.12/site-packages (from ipython>=7.23.1->ipykernel->jupyter) (2.15.1)
Collecting stack_data (from ipython>=7.23.1->ipykernel->jupyter)
  Downloading stack_data-0.6.3-py3-none-any.whl.metadata (18 kB)
Collecting wcwidth (from prompt_toolkit<3.1.0,>=3.0.41->ipython>=7.23.1->ipykernel->jupyter)
  Downloading wcwidth-0.2.13-py2.py3-none-any.whl.metadata (14 kB)
Collecting parso<0.9.0,>=0.8.4 (from jedi>=0.16->ipython>=7.23.1->ipykernel->jupyter)
  Downloading parso-0.8.4-py2.py3-none-any.whl.metadata (7.7 kB)
Requirement already satisfied: python-dateutil>=2.8.2 in /home/trian/.local/lib/python3.12/site-packages (from jupyter-client>=6.1.12->ipykernel->jupyter) (2.9.0.post0)
Requirement already satisfied: platformdirs>=2.5 in /home/trian/miniconda3/lib/python3.12/site-packages (from jupyter-core!=5.0.*,>=4.12->ipykernel->jupyter) (3.10.0)
Collecting ptyprocess>=0.5 (from pexpect>4.3->ipython>=7.23.1->ipykernel->jupyter)
  Downloading ptyprocess-0.7.0-py2.py3-none-any.whl.metadata (1.3 kB)
Requirement already satisfied: six>=1.5 in /home/trian/.local/lib/python3.12/site-packages (from python-dateutil>=2.8.2->jupyter-client>=6.1.12->ipykernel->jupyter) (1.17.0)
Collecting widgetsnbextension~=4.0.14 (from ipywidgets->jupyter)
  Downloading widgetsnbextension-4.0.14-py3-none-any.whl.metadata (1.6 kB)
Collecting jupyterlab_widgets~=3.0.15 (from ipywidgets->jupyter)
  Downloading jupyterlab_widgets-3.0.15-py3-none-any.whl.metadata (20 kB)
Collecting async-lru>=1.0.0 (from jupyterlab->jupyter)
  Downloading async_lru-2.0.5-py3-none-any.whl.metadata (4.5 kB)
Collecting httpx>=0.25.0 (from jupyterlab->jupyter)
  Downloading httpx-0.28.1-py3-none-any.whl.metadata (7.1 kB)
Collecting jinja2>=3.0.3 (from jupyterlab->jupyter)
  Downloading jinja2-3.1.6-py3-none-any.whl.metadata (2.9 kB)
Collecting jupyter-lsp>=2.0.0 (from jupyterlab->jupyter)
  Downloading jupyter_lsp-2.2.5-py3-none-any.whl.metadata (1.8 kB)
Collecting jupyter-server<3,>=2.4.0 (from jupyterlab->jupyter)
  Downloading jupyter_server-2.16.0-py3-none-any.whl.metadata (8.5 kB)
Collecting jupyterlab-server<3,>=2.27.1 (from jupyterlab->jupyter)
  Downloading jupyterlab_server-2.27.3-py3-none-any.whl.metadata (5.9 kB)
Collecting notebook-shim>=0.2 (from jupyterlab->jupyter)
  Downloading notebook_shim-0.2.4-py3-none-any.whl.metadata (4.0 kB)
Requirement already satisfied: setuptools>=41.1.0 in /home/trian/miniconda3/lib/python3.12/site-packages (from jupyterlab->jupyter) (75.8.0)
Collecting anyio>=3.1.0 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading anyio-4.9.0-py3-none-any.whl.metadata (4.7 kB)
Collecting argon2-cffi>=21.1 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading argon2_cffi-23.1.0-py3-none-any.whl.metadata (5.2 kB)
Collecting jupyter-events>=0.11.0 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading jupyter_events-0.12.0-py3-none-any.whl.metadata (5.8 kB)
Collecting jupyter-server-terminals>=0.4.4 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading jupyter_server_terminals-0.5.3-py3-none-any.whl.metadata (5.6 kB)
Collecting nbformat>=5.3.0 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading nbformat-5.10.4-py3-none-any.whl.metadata (3.6 kB)
Collecting overrides>=5.0 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading overrides-7.7.0-py3-none-any.whl.metadata (5.8 kB)
Collecting prometheus-client>=0.9 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading prometheus_client-0.22.0-py3-none-any.whl.metadata (14 kB)
Collecting send2trash>=1.8.2 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading Send2Trash-1.8.3-py3-none-any.whl.metadata (4.0 kB)
Collecting terminado>=0.8.3 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading terminado-0.18.1-py3-none-any.whl.metadata (5.8 kB)
Collecting websocket-client>=1.7 (from jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading websocket_client-1.8.0-py3-none-any.whl.metadata (8.0 kB)
Collecting babel>=2.10 (from jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading babel-2.17.0-py3-none-any.whl.metadata (2.0 kB)
Collecting json5>=0.9.0 (from jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading json5-0.12.0-py3-none-any.whl.metadata (36 kB)
Collecting jsonschema>=4.18.0 (from jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading jsonschema-4.23.0-py3-none-any.whl.metadata (7.9 kB)
Requirement already satisfied: requests>=2.31 in /home/trian/miniconda3/lib/python3.12/site-packages (from jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter) (2.32.3)
Requirement already satisfied: idna>=2.8 in /home/trian/miniconda3/lib/python3.12/site-packages (from anyio>=3.1.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter) (3.7)
Collecting sniffio>=1.1 (from anyio>=3.1.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading sniffio-1.3.1-py3-none-any.whl.metadata (3.9 kB)
Requirement already satisfied: typing_extensions>=4.5 in /home/trian/miniconda3/lib/python3.12/site-packages (from anyio>=3.1.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter) (4.12.2)
Collecting argon2-cffi-bindings (from argon2-cffi>=21.1->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading argon2_cffi_bindings-21.2.0-cp36-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (6.7 kB)
Requirement already satisfied: certifi in /home/trian/miniconda3/lib/python3.12/site-packages (from httpx>=0.25.0->jupyterlab->jupyter) (2025.1.31)
Collecting httpcore==1.* (from httpx>=0.25.0->jupyterlab->jupyter)
  Downloading httpcore-1.0.9-py3-none-any.whl.metadata (21 kB)
Collecting h11>=0.16 (from httpcore==1.*->httpx>=0.25.0->jupyterlab->jupyter)
  Downloading h11-0.16.0-py3-none-any.whl.metadata (8.3 kB)
Collecting MarkupSafe>=2.0 (from jinja2>=3.0.3->jupyterlab->jupyter)
  Downloading MarkupSafe-3.0.2-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (4.0 kB)
Collecting attrs>=22.2.0 (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading attrs-25.3.0-py3-none-any.whl.metadata (10 kB)
Collecting jsonschema-specifications>=2023.03.6 (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading jsonschema_specifications-2025.4.1-py3-none-any.whl.metadata (2.9 kB)
Collecting referencing>=0.28.4 (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading referencing-0.36.2-py3-none-any.whl.metadata (2.8 kB)
Collecting rpds-py>=0.7.1 (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter)
  Downloading rpds_py-0.25.1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (4.1 kB)
Collecting python-json-logger>=2.0.4 (from jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading python_json_logger-3.3.0-py3-none-any.whl.metadata (4.0 kB)
Collecting pyyaml>=5.3 (from jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading PyYAML-6.0.2-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (2.1 kB)
Collecting rfc3339-validator (from jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading rfc3339_validator-0.1.4-py2.py3-none-any.whl.metadata (1.5 kB)
Collecting rfc3986-validator>=0.1.1 (from jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading rfc3986_validator-0.1.1-py2.py3-none-any.whl.metadata (1.7 kB)
Collecting fqdn (from jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading fqdn-1.5.1-py3-none-any.whl.metadata (1.4 kB)
Collecting isoduration (from jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading isoduration-20.11.0-py3-none-any.whl.metadata (5.7 kB)
Requirement already satisfied: jsonpointer>1.13 in /home/trian/miniconda3/lib/python3.12/site-packages (from jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter) (2.1)
Collecting uri-template (from jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading uri_template-1.3.0-py3-none-any.whl.metadata (8.8 kB)
Collecting webcolors>=24.6.0 (from jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading webcolors-24.11.1-py3-none-any.whl.metadata (2.2 kB)
Collecting beautifulsoup4 (from nbconvert->jupyter)
  Downloading beautifulsoup4-4.13.4-py3-none-any.whl.metadata (3.8 kB)
Collecting bleach!=5.0.0 (from bleach[css]!=5.0.0->nbconvert->jupyter)
  Downloading bleach-6.2.0-py3-none-any.whl.metadata (30 kB)
Collecting defusedxml (from nbconvert->jupyter)
  Downloading defusedxml-0.7.1-py2.py3-none-any.whl.metadata (32 kB)
Collecting jupyterlab-pygments (from nbconvert->jupyter)
  Downloading jupyterlab_pygments-0.3.0-py3-none-any.whl.metadata (4.4 kB)
Collecting mistune<4,>=2.0.3 (from nbconvert->jupyter)
  Downloading mistune-3.1.3-py3-none-any.whl.metadata (1.8 kB)
Collecting nbclient>=0.5.0 (from nbconvert->jupyter)
  Downloading nbclient-0.10.2-py3-none-any.whl.metadata (8.3 kB)
Collecting pandocfilters>=1.4.1 (from nbconvert->jupyter)
  Downloading pandocfilters-1.5.1-py2.py3-none-any.whl.metadata (9.0 kB)
Collecting webencodings (from bleach!=5.0.0->bleach[css]!=5.0.0->nbconvert->jupyter)
  Downloading webencodings-0.5.1-py2.py3-none-any.whl.metadata (2.1 kB)
Collecting tinycss2<1.5,>=1.1.0 (from bleach[css]!=5.0.0->nbconvert->jupyter)
  Downloading tinycss2-1.4.0-py3-none-any.whl.metadata (3.0 kB)
Collecting fastjsonschema>=2.15 (from nbformat>=5.3.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading fastjsonschema-2.21.1-py3-none-any.whl.metadata (2.2 kB)
Requirement already satisfied: charset-normalizer<4,>=2 in /home/trian/miniconda3/lib/python3.12/site-packages (from requests>=2.31->jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter) (3.3.2)
Requirement already satisfied: urllib3<3,>=1.21.1 in /home/trian/miniconda3/lib/python3.12/site-packages (from requests>=2.31->jupyterlab-server<3,>=2.27.1->jupyterlab->jupyter) (2.3.0)
Requirement already satisfied: cffi>=1.0.1 in /home/trian/miniconda3/lib/python3.12/site-packages (from argon2-cffi-bindings->argon2-cffi>=21.1->jupyter-server<3,>=2.4.0->jupyterlab->jupyter) (1.17.1)
Requirement already satisfied: pycparser in /home/trian/miniconda3/lib/python3.12/site-packages (from cffi>=1.0.1->argon2-cffi-bindings->argon2-cffi>=21.1->jupyter-server<3,>=2.4.0->jupyterlab->jupyter) (2.21)
Collecting soupsieve>1.2 (from beautifulsoup4->nbconvert->jupyter)
  Downloading soupsieve-2.7-py3-none-any.whl.metadata (4.6 kB)
Collecting arrow>=0.15.0 (from isoduration->jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading arrow-1.3.0-py3-none-any.whl.metadata (7.5 kB)
Collecting types-python-dateutil>=2.8.10 (from arrow>=0.15.0->isoduration->jsonschema[format-nongpl]>=4.18.0->jupyter-events>=0.11.0->jupyter-server<3,>=2.4.0->jupyterlab->jupyter)
  Downloading types_python_dateutil-2.9.0.20250516-py3-none-any.whl.metadata (2.1 kB)
Collecting executing>=1.2.0 (from stack_data->ipython>=7.23.1->ipykernel->jupyter)
  Downloading executing-2.2.0-py2.py3-none-any.whl.metadata (8.9 kB)
Collecting asttokens>=2.1.0 (from stack_data->ipython>=7.23.1->ipykernel->jupyter)
  Downloading asttokens-3.0.0-py3-none-any.whl.metadata (4.7 kB)
Collecting pure-eval (from stack_data->ipython>=7.23.1->ipykernel->jupyter)
  Downloading pure_eval-0.2.3-py3-none-any.whl.metadata (6.3 kB)
Downloading jupyter-1.1.1-py2.py3-none-any.whl (2.7 kB)
Downloading plotly-6.1.1-py3-none-any.whl (16.1 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 16.1/16.1 MB 429.9 kB/s eta 0:00:00
Downloading narwhals-1.40.0-py3-none-any.whl (357 kB)
Downloading ipykernel-6.29.5-py3-none-any.whl (117 kB)
Downloading comm-0.2.2-py3-none-any.whl (7.2 kB)
Downloading debugpy-1.8.14-cp312-cp312-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl (4.2 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 4.2/4.2 MB 473.2 kB/s eta 0:00:00
Downloading ipython-9.2.0-py3-none-any.whl (604 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 604.3/604.3 kB 434.0 kB/s eta 0:00:00
Downloading prompt_toolkit-3.0.51-py3-none-any.whl (387 kB)
Downloading jedi-0.19.2-py2.py3-none-any.whl (1.6 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1.6/1.6 MB 348.1 kB/s eta 0:00:00
Downloading parso-0.8.4-py2.py3-none-any.whl (103 kB)
Downloading jupyter_client-8.6.3-py3-none-any.whl (106 kB)
Downloading jupyter_core-5.7.2-py3-none-any.whl (28 kB)
Downloading matplotlib_inline-0.1.7-py3-none-any.whl (9.9 kB)
Downloading pexpect-4.9.0-py2.py3-none-any.whl (63 kB)
Downloading ptyprocess-0.7.0-py2.py3-none-any.whl (13 kB)
Downloading pyzmq-26.4.0-cp312-cp312-manylinux_2_28_x86_64.whl (855 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 855.9/855.9 kB 307.6 kB/s eta 0:00:00
Downloading tornado-6.5.1-cp39-abi3-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl (443 kB)
Downloading traitlets-5.14.3-py3-none-any.whl (85 kB)
Downloading decorator-5.2.1-py3-none-any.whl (9.2 kB)
Downloading ipython_pygments_lexers-1.1.1-py3-none-any.whl (8.1 kB)
Downloading ipywidgets-8.1.7-py3-none-any.whl (139 kB)
Downloading jupyterlab_widgets-3.0.15-py3-none-any.whl (216 kB)
Downloading widgetsnbextension-4.0.14-py3-none-any.whl (2.2 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 2.2/2.2 MB 585.5 kB/s eta 0:00:00
Downloading jupyter_console-6.6.3-py3-none-any.whl (24 kB)
Downloading jupyterlab-4.4.2-py3-none-any.whl (12.3 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 12.3/12.3 MB 441.6 kB/s eta 0:00:00
Downloading jupyter_server-2.16.0-py3-none-any.whl (386 kB)
Downloading jupyterlab_server-2.27.3-py3-none-any.whl (59 kB)
Downloading anyio-4.9.0-py3-none-any.whl (100 kB)
Downloading argon2_cffi-23.1.0-py3-none-any.whl (15 kB)
Downloading async_lru-2.0.5-py3-none-any.whl (6.1 kB)
Downloading babel-2.17.0-py3-none-any.whl (10.2 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 10.2/10.2 MB 444.4 kB/s eta 0:00:00
Downloading httpx-0.28.1-py3-none-any.whl (73 kB)
Downloading httpcore-1.0.9-py3-none-any.whl (78 kB)
Downloading h11-0.16.0-py3-none-any.whl (37 kB)
Downloading jinja2-3.1.6-py3-none-any.whl (134 kB)
Downloading json5-0.12.0-py3-none-any.whl (36 kB)
Downloading jsonschema-4.23.0-py3-none-any.whl (88 kB)
Downloading attrs-25.3.0-py3-none-any.whl (63 kB)
Downloading jsonschema_specifications-2025.4.1-py3-none-any.whl (18 kB)
Downloading jupyter_events-0.12.0-py3-none-any.whl (19 kB)
Downloading jupyter_lsp-2.2.5-py3-none-any.whl (69 kB)
Downloading jupyter_server_terminals-0.5.3-py3-none-any.whl (13 kB)
Downloading MarkupSafe-3.0.2-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (23 kB)
Downloading nbconvert-7.16.6-py3-none-any.whl (258 kB)
Downloading mistune-3.1.3-py3-none-any.whl (53 kB)
Downloading bleach-6.2.0-py3-none-any.whl (163 kB)
Downloading tinycss2-1.4.0-py3-none-any.whl (26 kB)
Downloading nbclient-0.10.2-py3-none-any.whl (25 kB)
Downloading nbformat-5.10.4-py3-none-any.whl (78 kB)
Downloading fastjsonschema-2.21.1-py3-none-any.whl (23 kB)
Downloading notebook_shim-0.2.4-py3-none-any.whl (13 kB)
Downloading overrides-7.7.0-py3-none-any.whl (17 kB)
Downloading pandocfilters-1.5.1-py2.py3-none-any.whl (8.7 kB)
Downloading prometheus_client-0.22.0-py3-none-any.whl (62 kB)
Downloading python_json_logger-3.3.0-py3-none-any.whl (15 kB)
Downloading PyYAML-6.0.2-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (767 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 767.5/767.5 kB 658.9 kB/s eta 0:00:00
Downloading referencing-0.36.2-py3-none-any.whl (26 kB)
Downloading rfc3986_validator-0.1.1-py2.py3-none-any.whl (4.2 kB)
Downloading rpds_py-0.25.1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (390 kB)
Downloading Send2Trash-1.8.3-py3-none-any.whl (18 kB)
Downloading sniffio-1.3.1-py3-none-any.whl (10 kB)
Downloading terminado-0.18.1-py3-none-any.whl (14 kB)
Downloading webcolors-24.11.1-py3-none-any.whl (14 kB)
Downloading webencodings-0.5.1-py2.py3-none-any.whl (11 kB)
Downloading websocket_client-1.8.0-py3-none-any.whl (58 kB)
Downloading argon2_cffi_bindings-21.2.0-cp36-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (86 kB)
Downloading beautifulsoup4-4.13.4-py3-none-any.whl (187 kB)
Downloading soupsieve-2.7-py3-none-any.whl (36 kB)
Downloading defusedxml-0.7.1-py2.py3-none-any.whl (25 kB)
Downloading fqdn-1.5.1-py3-none-any.whl (9.1 kB)
Downloading isoduration-20.11.0-py3-none-any.whl (11 kB)
Downloading arrow-1.3.0-py3-none-any.whl (66 kB)
Downloading types_python_dateutil-2.9.0.20250516-py3-none-any.whl (14 kB)
Downloading jupyterlab_pygments-0.3.0-py3-none-any.whl (15 kB)
Downloading nest_asyncio-1.6.0-py3-none-any.whl (5.2 kB)
Downloading notebook-7.4.2-py3-none-any.whl (14.3 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 14.3/14.3 MB 333.7 kB/s eta 0:00:00
Downloading psutil-7.0.0-cp36-abi3-manylinux_2_12_x86_64.manylinux2010_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl (277 kB)
Downloading rfc3339_validator-0.1.4-py2.py3-none-any.whl (3.5 kB)
Downloading stack_data-0.6.3-py3-none-any.whl (24 kB)
Downloading asttokens-3.0.0-py3-none-any.whl (26 kB)
Downloading executing-2.2.0-py2.py3-none-any.whl (26 kB)
Downloading pure_eval-0.2.3-py3-none-any.whl (11 kB)
Downloading uri_template-1.3.0-py3-none-any.whl (11 kB)
Downloading wcwidth-0.2.13-py2.py3-none-any.whl (34 kB)
Installing collected packages: webencodings, wcwidth, pure-eval, ptyprocess, fastjsonschema, widgetsnbextension, websocket-client, webcolors, uri-template, types-python-dateutil, traitlets, tornado, tinycss2, soupsieve, sniffio, send2trash, rpds-py, rfc3986-validator, rfc3339-validator, pyzmq, pyyaml, python-json-logger, psutil, prompt_toolkit, prometheus-client, pexpect, parso, pandocfilters, overrides, nest-asyncio, narwhals, mistune, MarkupSafe, jupyterlab_widgets, jupyterlab-pygments, json5, ipython-pygments-lexers, h11, fqdn, executing, defusedxml, decorator, debugpy, bleach, babel, attrs, async-lru, asttokens, terminado, stack_data, referencing, plotly, matplotlib-inline, jupyter-core, jinja2, jedi, httpcore, comm, beautifulsoup4, arrow, argon2-cffi-bindings, anyio, jupyter-server-terminals, jupyter-client, jsonschema-specifications, isoduration, ipython, httpx, argon2-cffi, jsonschema, ipywidgets, ipykernel, nbformat, jupyter-console, nbclient, jupyter-events, nbconvert, jupyter-server, notebook-shim, jupyterlab-server, jupyter-lsp, jupyterlab, notebook, jupyter
Successfully installed MarkupSafe-3.0.2 anyio-4.9.0 argon2-cffi-23.1.0 argon2-cffi-bindings-21.2.0 arrow-1.3.0 asttokens-3.0.0 async-lru-2.0.5 attrs-25.3.0 babel-2.17.0 beautifulsoup4-4.13.4 bleach-6.2.0 comm-0.2.2 debugpy-1.8.14 decorator-5.2.1 defusedxml-0.7.1 executing-2.2.0 fastjsonschema-2.21.1 fqdn-1.5.1 h11-0.16.0 httpcore-1.0.9 httpx-0.28.1 ipykernel-6.29.5 ipython-9.2.0 ipython-pygments-lexers-1.1.1 ipywidgets-8.1.7 isoduration-20.11.0 jedi-0.19.2 jinja2-3.1.6 json5-0.12.0 jsonschema-4.23.0 jsonschema-specifications-2025.4.1 jupyter-1.1.1 jupyter-client-8.6.3 jupyter-console-6.6.3 jupyter-core-5.7.2 jupyter-events-0.12.0 jupyter-lsp-2.2.5 jupyter-server-2.16.0 jupyter-server-terminals-0.5.3 jupyterlab-4.4.2 jupyterlab-pygments-0.3.0 jupyterlab-server-2.27.3 jupyterlab_widgets-3.0.15 matplotlib-inline-0.1.7 mistune-3.1.3 narwhals-1.40.0 nbclient-0.10.2 nbconvert-7.16.6 nbformat-5.10.4 nest-asyncio-1.6.0 notebook-7.4.2 notebook-shim-0.2.4 overrides-7.7.0 pandocfilters-1.5.1 parso-0.8.4 pexpect-4.9.0 plotly-6.1.1 prometheus-client-0.22.0 prompt_toolkit-3.0.51 psutil-7.0.0 ptyprocess-0.7.0 pure-eval-0.2.3 python-json-logger-3.3.0 pyyaml-6.0.2 pyzmq-26.4.0 referencing-0.36.2 rfc3339-validator-0.1.4 rfc3986-validator-0.1.1 rpds-py-0.25.1 send2trash-1.8.3 sniffio-1.3.1 soupsieve-2.7 stack_data-0.6.3 terminado-0.18.1 tinycss2-1.4.0 tornado-6.5.1 traitlets-5.14.3 types-python-dateutil-2.9.0.20250516 uri-template-1.3.0 wcwidth-0.2.13 webcolors-24.11.1 webencodings-0.5.1 websocket-client-1.8.0 widgetsnbextension-4.0.14

# Generate dataset

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ python3 scripts/generate_sample_data.py
=== Fixed Single-Cell Data Generator ===
Generating synthetic data:
  Cells: 1,000
  Genes: 2,000
✅ Metadata saved to: data/synthetic/metadata.json

✅ Data generation completed!
   Output: data/raw/sparse_counts.tsv
   Entries: 278,568
   Genes: 2,000
   Cells: 1,000
   Total counts: 1,003,463
   Mean count: 3.60
   Sparsity: 0.861

Sample data preview:
 gene_idx  cell_idx     count
        0         3  2.000000
        0        10  3.285205
        0        15  6.087323
        0        27  4.000000
        0        45  2.169870
        0        50 14.311238
        0        54  3.000000
        0        63  1.000000
        0        66  1.000000
        0        70  3.614542

**data/synthetic/metadata.json**

{
  "total_entries": 278568,
  "unique_genes": 2000,
  "unique_cells": 1000,
  "total_counts": 1003462.6296762091,
  "mean_count": 3.6022178774166775,
  "median_count": 2.078111148222559,
  "max_count": 78.524232475965,
  "min_count": 1.0,
  "generation_parameters": {
    "n_cells": 1000,
    "n_genes": 2000,
    "n_cell_types": 5,
    "seed": 42
  },
  "sparsity": 0.860716
}

**data/raw/sparse_counts.tsv**

gene_idx	cell_idx	count
0	3	2.0
0	10	3.2852052488003753
0	15	6.087322986891296
0	27	4.0
0	45	2.1698702584588454
0	50	14.3112380407856
0	54	3.0
0	63	1.0
0	66	1.0
0	70	3.6145424047952215
0	71	1.0
0	83	1.0
0	85	2.4588632497636653
0	86	1.0
0	87	7.0
0	89	6.0
0	111	7.0
0	121	1.0
0	123	3.0
0	129	2.0
0	135	6.687081964827077
0	139	2.0
0	152	2.0
0	155	11.871807003547962
...

# Build Rust project

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ cargo build --release
   Compiling proc-macro2 v1.0.95
   Compiling libc v0.2.172
   Compiling memchr v2.7.4
   Compiling unicode-ident v1.0.18
   Compiling serde v1.0.219
   Compiling regex-syntax v0.8.5
   Compiling unicode-width v0.1.14
   Compiling vec_map v0.8.2
   Compiling ryu v1.0.20
   Compiling textwrap v0.11.0
   Compiling strsim v0.8.0
   Compiling bitflags v1.3.2
   Compiling ansi_term v0.12.1
   Compiling aho-corasick v1.1.3
   Compiling csv-core v0.1.12
   Compiling humantime v2.2.0
   Compiling itoa v1.0.15
   Compiling log v0.4.27
   Compiling termcolor v1.4.1
   Compiling regex-automata v0.4.9
   Compiling quote v1.0.40
   Compiling syn v2.0.101
   Compiling is-terminal v0.4.16
   Compiling atty v0.2.14
   Compiling clap v2.34.0
   Compiling regex v1.11.1
   Compiling serde_derive v1.0.219
   Compiling env_logger v0.10.2
   Compiling csv v1.3.1
   Compiling scrna-analyzer v1.0.0 (/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6)
    Finished `release` profile [optimized] target(s) in 54.29s

# Test with your data

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ ./target/release/scrna-analyzer \
    --input data/raw/sparse_counts.tsv \
    --output data/processed/cell_coords.tsv

**data/processed/cell_coords.tsv**

cell_id	pc1	pc2	pc3
0	-4.440650	-2.760000	-5.000000
1	-4.282455	-2.110000	-4.990000
2	-4.252528	-2.040000	-4.980000
3	-4.355582	-2.300000	-4.970000
4	-4.315052	-2.280000	-4.960000
5	-4.342045	-2.390000	-4.950000
6	-4.285194	-2.220000	-4.940000
7	-4.335899	-2.260000	-4.930000
8	-4.351956	-2.200000	-4.920000
9	-4.305184	-2.210000	-4.910000
10	-4.304658	-2.130000	-4.900000
11	-4.355916	-2.230000	-4.890000
12	-4.369058	-2.410000	-4.880000
13	-4.320234	-2.220000	-4.870000
14	-4.277820	-2.150000	-4.860000
15	-4.269160	-2.150000	-4.850000
16	-4.250170	-2.130000	-4.840000
17	-4.314302	-2.220000	-4.830000
18	-4.356547	-2.310000	-4.820000
19	-4.278464	-2.190000	-4.810000
20	-4.201784	-2.000000	-4.800000
...

# Run the Nextflow pipeline

(base) trian@triantoharyo:/mnt/c/Users/trian/BGVR/chapter_09/experiment_9_6$ nextflow run main.nf
Nextflow 25.04.2 is available - Please consider updating your version to it

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [magical_bartik] DSL2 - revision: 78435edf94


================================
Simple Single-Cell Pipeline
================================
Input: data/raw/sparse_counts.tsv
Output: results
================================

executor >  local (2)
[2f/dc4b39] VALIDATE_INPUT (1) [100%] 1 of 1 ✔
[08/b0620a] RUST_ANALYSIS (1)  [100%] 1 of 1 ✔
Pipeline completed: true

# results/coordinates/

**cell_coordinates.tsv**

cell_id	pc1	pc2	pc3
0	-4.440650	-2.760000	-5.000000
1	-4.282455	-2.110000	-4.990000
2	-4.252528	-2.040000	-4.980000
3	-4.355582	-2.300000	-4.970000
4	-4.315052	-2.280000	-4.960000
5	-4.342045	-2.390000	-4.950000
6	-4.285194	-2.220000	-4.940000
7	-4.335899	-2.260000	-4.930000
8	-4.351956	-2.200000	-4.920000
9	-4.305184	-2.210000	-4.910000
10	-4.304658	-2.130000	-4.900000
11	-4.355916	-2.230000	-4.890000
12	-4.369058	-2.410000	-4.880000
13	-4.320234	-2.220000	-4.870000
14	-4.277820	-2.150000	-4.860000
15	-4.269160	-2.150000	-4.850000
16	-4.250170	-2.130000	-4.840000
17	-4.314302	-2.220000	-4.830000
18	-4.356547	-2.310000	-4.820000
19	-4.278464	-2.190000	-4.810000
20	-4.201784	-2.000000	-4.800000
...

# results/pipeline_info/

**execution_trace.txt**

task_id	hash	native_id	name	status	exit	submit	duration	realtime	%cpu	peak_rss	peak_vmem	rchar	wchar
1	2f/dc4b39	63057	VALIDATE_INPUT (1)	COMPLETED	0	2025-05-25 19:15:12.689	2.2s	2s	83.0%	67.9 MB	430.3 MB	16.3 MB	4 MB
2	08/b0620a	63222	RUST_ANALYSIS (1)	COMPLETED	0	2025-05-25 19:15:14.989	515ms	359ms	26.5%	3.2 MB	4.6 MB	4.2 MB	32.8 KB

# results/validation/

**validated_matrix.tsv**

gene_idx	cell_idx	count
0	3	2.0
0	10	3.285205248800376
0	15	6.087322986891296
0	27	4.0
0	45	2.1698702584588454
0	50	14.3112380407856
0	54	3.0
0	63	1.0
0	66	1.0
0	70	3.6145424047952215
0	71	1.0
0	83	1.0
0	85	2.4588632497636653
0	86	1.0
0	87	7.0
0	89	6.0
0	111	7.0
0	121	1.0
0	123	3.0
0	129	2.0
0	135	6.687081964827077
0	139	2.0
0	152	2.0
0	155	11.871807003547962
...

**validation_report.txt**

Input Validation Report
======================
Original entries: 278,568
Non-zero entries: 278,568
Unique genes: 2,000
Unique cells: 1,000
Total counts: 1,003,463
```



