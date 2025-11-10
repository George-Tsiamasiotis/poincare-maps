use particle::Evolution;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use utils::{py_debug_impl, py_get_numpy1D, py_repr_impl};

#[derive(Clone)]
#[pyclass(name = "Evolution")]
pub struct PyEvolution(pub Evolution);

impl PyEvolution {
    /// Should only be created by PyParticle
    pub fn from_evolution(evolution: &Evolution) -> Self {
        Self(evolution.clone())
    }
}

#[pymethods]
impl PyEvolution {
    #[getter]
    pub fn get_duration(&self) -> String {
        format!("{:?}", self.0.duration)
    }

    #[getter]
    pub fn get_steps_taken(&self) -> usize {
        self.0.steps_taken()
    }

    #[getter]
    pub fn get_steps_stored(&self) -> usize {
        self.0.steps_stored()
    }
}

py_debug_impl!(PyEvolution);
py_repr_impl!(PyEvolution);
py_get_numpy1D!(PyEvolution, time);
py_get_numpy1D!(PyEvolution, theta);
py_get_numpy1D!(PyEvolution, psip);
py_get_numpy1D!(PyEvolution, rho);
py_get_numpy1D!(PyEvolution, zeta);
py_get_numpy1D!(PyEvolution, psi);
py_get_numpy1D!(PyEvolution, ptheta);
py_get_numpy1D!(PyEvolution, pzeta);
py_get_numpy1D!(PyEvolution, energy);
