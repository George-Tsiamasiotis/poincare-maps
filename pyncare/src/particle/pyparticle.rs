use crate::{PyBfield, PyCurrents, PyPerturbation, PyQfactor};
use crate::{PyEvolution, PyMappingParameters, PyParticleError};
use particle::{InitialConditions, Particle};
use utils::{py_debug_impl, py_repr_impl};

use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "InitialConditions")]
pub struct PyInitialConditions(pub InitialConditions);

#[pymethods]
impl PyInitialConditions {
    /// Creates a new PyInitialConditions wrapper.
    #[new]
    #[rustfmt::skip]
    pub fn new(time0: f64, theta0: f64, psip0: f64, rho0: f64, zeta0: f64, mu: f64) -> Self {
        Self (InitialConditions {time0, theta0, psip0, rho0, zeta0, mu})
    }
}

py_debug_impl!(PyInitialConditions);
py_repr_impl!(PyInitialConditions);

#[pyclass(name = "Particle")]
pub struct PyParticle(pub Particle);

#[pymethods]
impl PyParticle {
    /// Creates a new PyParticle wrapper.
    #[new]
    pub fn new(initial: &PyInitialConditions) -> Self {
        Self(Particle::new(&initial.0))
    }

    pub fn integrate<'py>(
        &mut self,
        qfactor: &PyQfactor,
        bfield: &PyBfield,
        currents: &PyCurrents,
        perturbation: &PyPerturbation,
        t_eval: Bound<'py, PyList>,
    ) -> Result<(), PyParticleError> {
        // This is necessary since t_eval is a sequence
        let t_eval_vec: Vec<f64> = t_eval.iter().map(|t| t.extract().unwrap()).collect();
        match t_eval_vec.len() {
            2 => (),
            _ => panic!("`t_eval` must be of the form (t0, tf)"),
        };

        let t_eval: (f64, f64) = (t_eval_vec[0], t_eval_vec[1]);
        self.0
            .integrate(&qfactor.0, &bfield.0, &currents.0, &perturbation.0, t_eval)?;
        Ok(())
    }

    pub fn map(
        &mut self,
        qfactor: &PyQfactor,
        bfield: &PyBfield,
        currents: &PyCurrents,
        perturbation: &PyPerturbation,
        mapping: &PyMappingParameters,
    ) -> Result<(), PyParticleError> {
        Ok(self.0.map(
            &qfactor.0,
            &bfield.0,
            &currents.0,
            &perturbation.0,
            &mapping.0,
        )?)
    }

    #[getter]
    pub fn get_evolution(&self) -> PyEvolution {
        PyEvolution::from_evolution(&self.0.evolution)
    }

    #[getter]
    pub fn get_status(&self) -> String {
        format!("{:?}", self.0.status)
    }
}

py_debug_impl!(PyParticle);
py_repr_impl!(PyParticle);
