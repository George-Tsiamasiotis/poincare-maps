mod error;
mod particle;
mod pylibrium;

use pyo3::prelude::*;

pub use particle::*;
pub use pylibrium::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Pylibrium
    m.add_class::<PyQfactor>()?;
    m.add_class::<PyCurrent>()?;
    m.add_class::<PyBfield>()?;
    m.add_class::<PyHarmonic>()?;
    m.add_class::<PyPerturbation>()?;
    // Particle
    m.add_class::<PyInitialConditions>()?;
    m.add_class::<PyParticle>()?;
    m.add_class::<PyEvolution>()?;
    Ok(())
}
