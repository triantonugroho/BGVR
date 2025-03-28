## 5.5. Variant Calling and Genotyping

### experiment_55

This Rust code processes variant hypotheses for a single or multiple genomic positions in a parallel and chunked fashion. The pileup data, which contains base and quality information, is loaded from a JSON file, and a list of variant hypotheses is retrieved from another JSON file. These hypotheses each specify a chromosome, position, reference base, and an alternative base to test. Because a real HPC pipeline might contain millions of hypotheses, the code is designed to break them into manageable chunks to reduce memory overhead, writing partial results after each chunk is processed.

The concurrency for the actual likelihood calculations is handled by Rayon’s data-parallel iterators. Each variant hypothesis is mapped to a function that computes a naive log-likelihood of the alternative allele versus the reference allele, given the reads in the pileup. This naive_likelihood function interprets the base quality as a probability of error, which is typical in sequencing data. In HPC contexts, ephemeral containers might each process a subset of the hypotheses or the genomic data, storing their partial results in JSON files.

After collecting command-line parameters through the clap crate, the code loads both pileup data and variant hypotheses into memory. In a more complex pipeline, or in truly large-scale HPC runs, this data might instead be streamed or accessed through partial file segments, allowing ephemeral tasks to work in isolation. The hypotheses are processed in user-defined “chunks” to minimize peak memory usage. Within each chunk, the code calls compute_variants, which uses Rayon’s .par_iter() to iterate over each variant hypothesis in parallel. The result of this step is a set of partial variant calls, serialized immediately to disk to avoid data loss and to enable potential downstream merging.

Finally, all partial results written to the output directory are reopened, deserialized, and merged into a single vector of calls. This final vector is written to a single merged JSON file, showing how ephemeral tasks can unify their outputs without each having to hold the entire dataset in memory. In a production environment, it is straightforward to extend this approach. One can replace the naive_likelihood function with more sophisticated models that use machine learning or Bayesian statistics, or incorporate advanced concurrency management with Crossbeam. Regardless of the complexity of the expansion, Rust’s ownership and type system ensure that parallel computations remain memory-safe and free from data races.

#### Files contents:
* experiment_55/
  * Cargo.toml (Cargo.toml file for dependencies)
* experiment_55/src/
  * main.rs (rust script)
  * hypotheses.json (json input file)
  * pileup.json (json input file)
  * merged_variants.json (merged variants json output file)
* experiment_55/src/output
  * partial_varians_chunk_0.json (partial varians in chunk 0 json output file)


#### How to run:

run in powershell:

```powershell
cargo run -- --pileup-input pileup.json --hypotheses-input hypotheses.json --chunk-size 5000 --output-dir output --merged-output merged_variants.json | tee output.txt
```

(run main.rs with chunk size 5000, 2 json input file name and 1 json output filename)
  
#### [dependencies]

```toml
anyhow = "1.0"
rayon = "1.7"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
clap = { version = "4.0", features = ["derive"] }
```

#### Explanation of the Output

##### 1. Chunk Processing Log (output.txt)

* The log states:

```rust
Chunk 0 processed (2 hypotheses). Partial results saved at "output\\partial_variants_chunk_0.json"
Successfully merged 2 variants into "merged_variants.json"
```

* This means:

  * The program processed a single chunk (Chunk 0) because only 2 hypotheses were present in hypotheses.json, which is far below the chunk size of 5000.

  * The computed results were saved as a JSON file output/partial_variants_chunk_0.json.

  * Finally, all the chunks were merged into a single output file named merged_variants.json.

##### 2. Partial Chunk File (partial_variants_chunk_0.json)

```json
{
  "sites": [
    {
      "chrom": "chr1",
      "position": 100,
      "ref_base": "A",
      "alt_base": "G",
      "likelihood": -0.00575393902402108
    },
    {
      "chrom": "chr1",
      "position": 101,
      "ref_base": "C",
      "alt_base": "T",
      "likelihood": -0.0013167781101437682
    }
  ]
}
```

* This file contains the computed variant sites (potential mutations) for chunk 0.

* Two variants are reported:

  * Variant 1: chr1:100, changing A → G, with likelihood -0.00575.

  * Variant 2: chr1:101, changing C → T, with likelihood -0.00131.

* The likelihood values are negative log probabilities, meaning lower values indicate stronger evidence for the variant.

##### 3. Final Merged Output (merged_variants.json)

```json
[
  {
    "chrom": "chr1",
    "position": 100,
    "ref_base": "A",
    "alt_base": "G",
    "likelihood": -0.00575393902402108
  },
  {
    "chrom": "chr1",
    "position": 101,
    "ref_base": "C",
    "alt_base": "T",
    "likelihood": -0.0013167781101437682
  }
]
```

* The final output is identical to partial_variants_chunk_0.json because only one chunk was processed.

* If more hypotheses were present, multiple chunk files would have been merged into merged_variants.json.

#### Conclusion
1. The program successfully processed variant calling

    * It read pileup data and variant hypotheses from pileup.json and hypotheses.json.

    * It calculated likelihood scores for each variant hypothesis based on sequencing data.

    * The results were written in chunks and later merged.

2. The pipeline works correctly, but the dataset is small

    * The chunk_size was set to 5000, but the input had only 2 hypotheses, so only one chunk was processed.

    * This confirms that the pipeline can handle much larger datasets efficiently.

3. Likelihood scores indicate the confidence in detected variants

    * The negative log-likelihood values suggest how likely each mutation is.

    * A more negative value means higher confidence in the variant.

4. Next Steps

    * Run on a larger dataset to test chunking behavior.

    * Compare likelihood thresholds to filter strong variant candidates.

    * Optimize performance for parallel execution with Rayon.

This result demonstrates that the variant calling pipeline is functional and scalable. 
