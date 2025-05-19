## 2.5. Introduction Machine Learning in Bioinformatics

### experiment_2_5

Below is an illustrative Rust program that demonstrates how one might employ a deep neural network, via the tch-rs crate (Rust bindings for PyTorch), to classify genomic or transcriptomic data. In this hypothetical scenario, we assume we have a CSV file containing gene expression values along with sample labels indicating, for example, the presence or absence of a particular disease phenotype. Our goal is to train a multi-layer perceptron (MLP) to distinguish positive from negative cases based on their expression profiles.

In this example, “gene_expression.csv” contains rows of numerical values—one row per sample—and a final column indicating the label (e.g., 0 or 1). The program reads the CSV, builds a simple MLP, and iterates through mini-batches to train a binary classifier. While this code focuses on a single-file dataset, in real-world scenarios you may have gigabytes of multi-omics data split across multiple files, and you might integrate concurrency or distributed computing patterns for scalability. Nevertheless, it illustrates the essential building blocks of using tch-rs in bioinformatics.

#### Project Structure:

```plaintext
experiment_2_5/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    ├── Python Code for Synthesize gene_expression.csv.ipynb   # Python notebook for data synthesis
    ├── gene_expression.csv        # Gene expression CSV file (Google Drive link)
    └── output.txt                 # Output file
```
  
#### Cargo.toml

```toml
[package]
name = "experiment_2_5"
version = "0.1.0"
edition = "2021"

[dependencies]
tch = "0.19.0"
```

#### How to run

```powershell
cargo run | tee output.txt
```
(run main.rs and save the output in output.txt)
  

#### Explanation of the Output

The output represents the training process of a neural network model using gene expression data from a CSV file. The model is trained for 10 epochs, and the validation loss is recorded at each epoch. Below is a breakdown of the key aspects of the output:

### 1. Training Process

* The model is a fully connected neural network (FCNN) with three layers:
  * Input layer: 10,000 features (gene expression values).
  * Hidden layers: 128 and 64 neurons with ReLU activation.
  * Output layer: 1 neuron (for binary classification).
* The loss function used is Binary Cross-Entropy with Logits Loss (BCEWithLogitsLoss).
*The optimizer is Adam with a learning rate of 1e-3.
* Training is performed in mini-batches of 64 samples.

##### 2. Validation Loss Progression

Epoch 1, Validation loss: 71.8563

Epoch 2, Validation loss: 59.6753

Epoch 3, Validation loss: 53.6821

Epoch 4, Validation loss: 51.4969

Epoch 5, Validation loss: 46.6416

Epoch 6, Validation loss: 47.4496

Epoch 7, Validation loss: 38.8133

Epoch 8, Validation loss: 36.7464

Epoch 9, Validation loss: 31.4197

Epoch 10, Validation loss: 25.9429

* The validation loss starts at 71.8563 in Epoch 1 and gradually decreases, reaching 25.9429 in Epoch 10.
* A decreasing loss indicates that the model is learning and improving in predicting the labels.
* There is a slight increase in loss at Epoch 6 (47.4496), which could be due to:
  * Noisy data or an unstable learning rate.
  * Overfitting to the training data.
  * A difficult batch of samples.
* After Epoch 6, the loss continues to decrease, suggesting that the model is still improving.
  
##### 3. Key Insights

* Effective Learning: The steady decline in validation loss suggests that the model is effectively learning from the dataset.
* Potential for Further Training: Since the loss is still decreasing, additional epochs could further improve performance.
* Possible Overfitting Check: If the loss plateaus or increases after many epochs, regularization techniques like dropout or early stopping may be needed.
* Evaluation Required: Final model performance should be validated using a test dataset to ensure it generalizes well.
  
#### Conclusion
The neural network successfully learns from the gene expression dataset, significantly reducing validation loss over 10 epochs. This suggests that the model is effectively capturing patterns in the data and improving its predictions.




