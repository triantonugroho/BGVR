use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

/// Fungsi Reverse Complement DNA
#[pyfunction]
fn reverse_complement(seq: &str) -> PyResult<String> {
    let rev_comp: String = seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T', 'T' => 'A',
            'C' => 'G', 'G' => 'C',
            _ => c,
        })
        .collect();
    Ok(rev_comp)
}

/// Menghitung GC Content dalam persen
#[pyfunction]
fn gc_content(seq: &str) -> PyResult<f64> {
    if seq.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("Input sequence cannot be empty"));
    }
    let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count();
    let gc_percentage = (gc_count as f64 / seq.len() as f64) * 100.0;
    Ok(gc_percentage)
}

/// Modul Python
#[pymodule]
fn pyo3_biopython_example(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(gc_content, m)?)?;
    Ok(())
}
