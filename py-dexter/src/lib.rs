mod args;
mod equilibrium;
mod error;
mod macros;

pub use args::*;
pub use dexter::dexter_equilibrium::*;
pub use equilibrium::*;
pub use error::*;

pub type Result<T> = std::result::Result<T, DexterError>;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<equilibrium::PyLastClosedFluxSurface>()?;
    m.add_class::<equilibrium::PyQfactor>()?;
    m.add_class::<equilibrium::PyCurrent>()?;
    Ok(())
}
