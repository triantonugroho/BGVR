## 1.4. AI/ML Implementation with Rust Crates

### experiment_1.4

Below is a sample Rust code snippet demonstrating how to do K-means clustering with linfa and then, in a separate pipeline, train a neural network with tch-rs. Each step includes commentary on how to handle industrial-scale issues like concurrency and data loading.

#### Project Structure:

```plaintext
experiment_1.4/
├── Cargo.toml                     # Rust project configuration and dependencies
└── src/
    ├── main.rs                    # Main Rust script containing program logic
    └── output.txt                 # Output file
```

#### How to run:

run in powershell:

```powershell
cargo run | tee output.txt
```

(run main.rs and save the output in output.txt)
  
#### Cargo.toml

```toml
[package]
name = "experiment_1.4"
version = "0.1.0"
edition = "2021"

[dependencies]
linfa = "0.7.1"
linfa-clustering = "0.7.1"
tch = "0.19.0"
anyhow = "1.0"
thiserror = "2.0.11"
ndarray = { version = "0.15", features = ["approx"] }
```

#### Explanation of the Output

This Rust program demonstrates two machine learning techniques using the linfa and tch (Torch for Rust) libraries:
1. K-Means Clustering – Unsupervised learning technique for clustering data.
2. Neural Network Training – A simple feedforward neural network for classification.

##### Step-by-Step Execution

###### 1. K-Means Clustering
The function run_kmeans_example() performs the following steps:

####### Step 1: Define Data

```rust
let data: Array2<f64> = ndarray::arr2(&[
    [2.1, 3.4],
    [2.0, 3.2],
    [8.9, 9.1],
    [9.0, 8.7],
]);
```

* The dataset consists of 4 points in 2D space.
* The goal is to cluster these points into 2 groups.

####### Step 2: Train K-Means Model

```rust
let model = KMeans::params(2)
    .max_n_iterations(100)
    .fit(&dataset)?;
```

* The model is initialized with 2 clusters (K = 2).
* It iterates 100 times to find the best cluster assignments.

####### Step 3: Predict Clusters

```rust
let labels = model.predict(&dataset);
println!("K-means cluster assignments: {:?}", labels);
```

* The model assigns each point to a cluster.
  
####### Output:

```sh
K-means cluster assignments: [0, 0, 1, 1], shape=[4], strides=[1], layout=CFcf (0xf), const ndim=1
```

* The first two points ([2.1, 3.4] and [2.0, 3.2]) are in Cluster 0.
* The last two points ([8.9, 9.1] and [9.0, 8.7]) are in Cluster 1.
* This makes sense because the two groups of points are clearly separate in the 2D space.

###### 2. Neural Network Training
The function run_nn_example() implements a simple feedforward neural network.

####### Step 1: Define Neural Network
```rust
let net = {
    let root = &vs.root();
    nn::seq()
        .add(nn::linear(root, 10, 5, Default::default())) // Input: 10 features → 5 hidden units
        .add_fn(|xs| xs.relu())                          // ReLU activation
        .add(nn::linear(root, 5, 2, Default::default()))  // Hidden: 5 → Output: 2 classes
};
```

* Input layer: 10 features per sample.
* Hidden layer: 5 neurons, ReLU activation.
* Output layer: 2 neurons (for binary classification).

####### Step 2: Generate Training Data

```rust
let x_train = Tensor::rand(&[100, 10], (tch::Kind::Float, device));
let y_train = Tensor::randint(0, 2, (tch::Kind::Int64, device));
```

* x_train: 100 samples, each with 10 random features.
* y_train: 100 random labels (0 or 1).

####### Step 3: Train the Model

```rust
for epoch in 1..101 {
    let logits = net.forward(&x_train);
    let loss = logits.cross_entropy_for_logits(&y_train);
    opt.backward_step(&loss);

    if epoch % 20 == 0 {
        println!("Epoch: {}, loss: {:?}", epoch, f64::from(loss.double_value(&[0])));
    }
}
```

* Runs for 100 epochs.
* Uses cross-entropy loss for classification.
* Prints the loss every 20 epochs.

#### Conclusion

##### 1. K-Means Clustering

* Successfully grouped 4 data points into 2 clusters.
* The output confirms that the algorithm correctly detected two distinct groups.

##### 2. Neural Network

* Defined and trained a simple feedforward neural network.
* Uses SGD optimizer and cross-entropy loss.
* The model would improve over training epochs.

This program demonstrates both unsupervised (clustering) and supervised (neural network) learning using Rust!
