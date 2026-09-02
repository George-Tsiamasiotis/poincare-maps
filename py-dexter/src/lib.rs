mod args;
mod equilibrium;
mod error;
mod macros;
mod simulate;

pub use args::*;
pub use equilibrium::*;
pub use error::*;

pub type Result<T> = std::result::Result<T, DexterError>;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<equilibrium::PyLastClosedFluxSurface>()?;
    m.add_class::<equilibrium::PyGeometry>()?;
    m.add_class::<equilibrium::PyQfactor>()?;
    m.add_class::<equilibrium::PyCurrent>()?;
    m.add_class::<equilibrium::PyBfield>()?;
    m.add_class::<equilibrium::PyMode>()?;
    m.add_class::<equilibrium::PyPerturbation>()?;
    m.add_class::<simulate::PyInitialFlux>()?;
    m.add_class::<simulate::PyInitialConditions>()?;
    Ok(())
}
