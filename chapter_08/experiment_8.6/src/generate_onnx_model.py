# This is a simplified representation of an ONNX model.
# In a real scenario, you would generate a proper .onnx file using frameworks like PyTorch or TensorFlow.
# For our demonstration purposes, I'll create a script to generate a very simple ONNX model.

import onnx
from onnx import TensorProto
from onnx.helper import make_model, make_node, make_graph, make_tensor_value_info, make_tensor
import numpy as np

# Create a simple model that adds the features and applies a sigmoid function
# Input: [ref_length, alt_length, node_degree, centrality, complexity]
# Output: probability score

# Define the graph
X = make_tensor_value_info('input', TensorProto.FLOAT, [None, 5])  # Batch size x 5 features
Y = make_tensor_value_info('output', TensorProto.FLOAT, [None, 1])  # Batch size x 1 score

# Define the weights (these would be learned weights in a real model)
weight_tensor = make_tensor('weights', TensorProto.FLOAT, [5, 1], 
                           np.array([0.2, -0.1, 0.15, 0.3, 0.25], dtype=np.float32).reshape(5, 1).tobytes(), 
                           raw=True)
bias_tensor = make_tensor('bias', TensorProto.FLOAT, [1], 
                         np.array([0.1], dtype=np.float32).tobytes(), 
                         raw=True)

# Define nodes
matmul_node = make_node('MatMul', ['input', 'weights'], ['matmul_output'])
add_node = make_node('Add', ['matmul_output', 'bias'], ['add_output'])
sigmoid_node = make_node('Sigmoid', ['add_output'], ['output'])

# Create the graph
graph = make_graph(
    [matmul_node, add_node, sigmoid_node],
    'variant-scoring-model',
    [X],
    [Y],
    [weight_tensor, bias_tensor]
)

# Create the model
model = make_model(graph, producer_name='variant-scorer-sample')

# Save the model
onnx.save(model, 'variant_model.onnx')

print("Sample ONNX model saved to variant_model.onnx")