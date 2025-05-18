use pyo3::{prelude::*, exceptions::PyValueError, wrap_pyfunction};
use tch::{nn, nn::Module, Tensor, Device};
use std::convert::TryInto;

#[pyfunction]
fn run_linear_forward_py(input_data: Vec<f32>) -> PyResult<Vec<f32>> {
    // Pastikan input memiliki dua elemen
    if input_data.len() != 2 {
        return Err(PyErr::new::<PyValueError, _>(
            "Input data must contain exactly two float values.",
        ));
    }

    // Gunakan CPU sebagai perangkat eksekusi
    let device = Device::Cpu;
    let vs = nn::VarStore::new(device);
    let root = vs.root();

    // Buat layer linear dengan 2 input dan 1 output
    let linear = nn::linear(&root, 2, 1, Default::default());

    // Konversi input menjadi tensor 1x2
    let input_tensor = Tensor::from_slice(&input_data)
        .view([1, 2])
        .to_device(device);

    // Jalankan forward pass
    let output = linear.forward(&input_tensor);

    // Konversi output tensor menjadi Vec<f32>
    let out_vec: Vec<f32> = output.try_into().map_err(|e| {
        PyErr::new::<PyValueError, _>(format!("Tensor conversion error: {}", e))
    })?;

    Ok(out_vec)
}

#[pymodule]
fn pyo3_tch_example(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Registrasi fungsi agar bisa diakses dari Python
    m.add_function(wrap_pyfunction!(run_linear_forward_py, py)?)?;
    Ok(())
}
