mod equilibrium;
mod error;
mod initial;
mod point;
mod rk45;
mod state;
mod system;

pub use equilibrium::Bfield;
pub use equilibrium::Current;
pub use equilibrium::Qfactor;

pub use error::MapError;
pub use initial::InitialConditions;
pub use system::System;

pub(crate) use point::Point;
pub(crate) use rk45::Rk45State;
pub(crate) use state::State;

pub type Result<T> = std::result::Result<T, MapError>;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use pyo3::prelude::*;

#[pymodule]
fn poincare_maps(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InitialConditions>()?;
    m.add_class::<Bfield>()?;
    m.add_class::<Qfactor>()?;
    m.add_class::<Current>()?;
    Ok(())
}
