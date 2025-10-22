#![feature(coverage_attribute)]

mod consts;
mod equilibrium;
mod error;
mod evolution;
mod particle;
mod poincare;
mod point;
mod solver;
mod state;
mod utils;

pub use evolution::Evolution;
pub use point::Point;

pub use equilibrium::Bfield;
pub use equilibrium::Current;
pub use equilibrium::Harmonic;
pub use equilibrium::Perturbation;
pub use equilibrium::Qfactor;

pub use error::MapError;
pub use particle::Particle;
pub use poincare::Poincare;
pub use state::State;

pub use consts::*;
pub use utils::*;

pub type Result<T> = std::result::Result<T, MapError>;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use pyo3::prelude::*;

#[pymodule]
#[coverage(off)]
fn poincare_maps(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Bfield>()?;
    m.add_class::<Qfactor>()?;
    m.add_class::<Current>()?;
    m.add_class::<Harmonic>()?;
    m.add_class::<Perturbation>()?;
    m.add_class::<Particle>()?;
    m.add_class::<Poincare>()?;
    m.add_class::<PoincareParameters>()?;
    Ok(())
}
