mod args;
mod error;
mod machine;
mod macros;
mod simulate;

pub use args::*;
pub use error::*;
pub use machine::*;
pub use simulate::*;

pub type Result<T> = std::result::Result<T, DexterError>;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<machine::PyLastClosedFluxSurface>()?;
    m.add_class::<machine::PyGeometry>()?;
    m.add_class::<machine::PyQfactor>()?;
    m.add_class::<machine::PyCurrent>()?;
    m.add_class::<machine::PyBfield>()?;
    m.add_class::<machine::PyMode>()?;
    m.add_class::<machine::PyPerturbation>()?;
    m.add_class::<simulate::PyInitialFlux>()?;
    m.add_class::<simulate::PyInitialConditions>()?;
    Ok(())
}
