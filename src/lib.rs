mod equilibrium;
mod error;
mod initial;
mod particle;
mod solver;
mod state;

pub use equilibrium::Bfield;
pub use equilibrium::Current;
pub use equilibrium::Qfactor;

pub use error::MapError;
pub use initial::InitialConditions;
pub use particle::Particle;

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
    m.add_class::<Particle>()?;
    Ok(())
}
