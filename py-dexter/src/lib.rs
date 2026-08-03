mod equilibrium;
mod macros;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<equilibrium::PyLastClosedFluxSurface>()?;
    Ok(())
}
