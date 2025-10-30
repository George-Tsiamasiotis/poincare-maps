use particle::Particle;

use pyo3::prelude::*;
use pyo3::types::PyList;

use crate::error::PyParticleError;
use crate::{repr_impl, PyBfield, PyCurrent, PyPerturbation, PyQfactor};
use crate::{PyEvolution, PyInitialConditions};

#[pyclass]
#[pyo3(name = "Particle")]
pub struct PyParticle {
    pub particle: Particle,
}

#[pymethods]
impl PyParticle {
    #[new]
    pub fn new(py_initial: &PyInitialConditions) -> Self {
        let particle = Particle::new(&py_initial.initial);
        Self { particle }
    }

    pub fn integrate<'py>(
        &mut self,
        qfactor: &PyQfactor,
        bfield: &PyBfield,
        current: &PyCurrent,
        per: &PyPerturbation,
        t_eval: Bound<'py, PyList>,
    ) -> Result<(), PyParticleError> {
        // This is necessary since t_eval is a sequence
        let t_eval_vec: Vec<f64> = t_eval.iter().map(|t| t.extract().unwrap()).collect();
        match t_eval_vec.len() {
            2 => (),
            _ => panic!("`t_eval`, must be of the form (t0, tf)"),
        };

        let t_eval: (f64, f64) = (t_eval_vec[0], t_eval_vec[1]);
        self.particle.integrate(
            &qfactor.qfactor,
            &bfield.bfield,
            &current.current,
            &per.perturbation,
            t_eval,
        )?;
        Ok(())
    }

    // TODO: add mapping and test it

    #[getter]
    pub fn get_evolution(&self) -> PyEvolution {
        PyEvolution::from_evolution(&self.particle.evolution)
    }

    #[getter]
    pub fn get_status(&self) -> String {
        format!("{:?}", self.particle.status)
    }
}

repr_impl!(PyParticle);

impl std::fmt::Debug for PyParticle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.particle.fmt(f)
    }
}
