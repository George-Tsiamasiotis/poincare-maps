mod error;
mod pylibrium;

use pyo3::prelude::*;

pub use pylibrium::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyQfactor>()?;
    m.add_class::<PyCurrent>()?;
    Ok(())
}
