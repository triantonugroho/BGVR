use linfa::prelude::*;
use linfa_clustering::KMeans;
use ndarray::Array2;
use tch::{nn, nn::Module, nn::OptimizerConfig, Device, Tensor};

fn run_kmeans_example() -> Result<(), Box<dyn std::error::Error>> {
    // Define the 2D data (e.g., 2 features per sample)
    let data: Array2<f64> = ndarray::arr2(&[
        [2.1, 3.4],
        [2.0, 3.2],
        [8.9, 9.1],
        [9.0, 8.7],
    ]);

    // KMeans does not need target labels, so we only provide the data
    let dataset = linfa::Dataset::from(data);
    
    // Fit a KMeans model with 2 clusters
    let model = KMeans::params(2)
        .max_n_iterations(100)
        .fit(&dataset)?;

    // Get cluster assignments for each point
    let labels = model.predict(&dataset);
    println!("K-means cluster assignments: {:?}", labels);

    Ok(())
}

fn run_nn_example() -> Result<(), Box<dyn std::error::Error>> {
    let device = Device::cuda_if_available(); // Use CUDA if available, otherwise fall back to CPU
    let vs = nn::VarStore::new(device);

    // Define a simple neural network
    let net = {
        let root = &vs.root();
        nn::seq()
            .add(nn::linear(root, 10, 5, Default::default())) // 10 input features, 5 hidden units
            .add_fn(|xs| xs.relu())                          // ReLU activation
            .add(nn::linear(root, 5, 2, Default::default()))  // 5 hidden units, 2 output classes
    };

    // Create an SGD optimizer
    let mut opt = nn::Sgd::default().build(&vs, 1e-3)?;

    // Generate training data: tensors of size [100, 10] for features and [100] for labels
    let x_train = Tensor::rand(&[100, 10], (tch::Kind::Float, device)); // 100 samples with 10 features each
    let y_train = Tensor::randint(0, 2, (tch::Kind::Int64, device)); // Binary labels (0 or 1) for 100 samples
        
    // Training loop (100 epochs)
    for epoch in 1..101 {
        let logits = net.forward(&x_train); // Forward pass
        let loss = logits.cross_entropy_for_logits(&y_train); // Cross-entropy loss for classification
        opt.backward_step(&loss); // Backward pass and optimization

        // Print loss every 20 epochs
        if epoch % 20 == 0 {
            println!("Epoch: {}, loss: {:?}", epoch, f64::from(loss.double_value(&[0])));
        }
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Run KMeans clustering example
    run_kmeans_example()?;

    // Run neural network training example
    run_nn_example()?;

    Ok(())
}
