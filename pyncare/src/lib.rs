mod particle;
mod poincare;
mod pyerrors;
mod pylibrium;

use pyo3::prelude::*;

pub use particle::*;
pub use poincare::*;
pub use pyerrors::*;
pub use pylibrium::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Pylibrium
    m.add_class::<PyQfactor>()?;
    m.add_class::<PyCurrents>()?;
    m.add_class::<PyBfield>()?;
    m.add_class::<PyHarmonic>()?;
    m.add_class::<PyPerturbation>()?;
    // Particle
    m.add_class::<PyInitialConditions>()?;
    m.add_class::<PyParticle>()?;
    m.add_class::<PyEvolution>()?;
    m.add_class::<PyMapping>()?;
    // Poincare
    m.add_class::<PyPoincareInit>()?;
    m.add_class::<PyPoincare>()?;
    Ok(())
}
