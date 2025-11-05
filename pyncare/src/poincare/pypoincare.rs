use poincare::Poincare;
use utils::to_numpy2D_impl;

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use crate::{PyBfield, PyCurrent, PyMapping, PyPerturbation, PyQfactor};
use crate::{PyPoincareError, PyPoincareInit};

#[pyclass(name = "Poincare")]
pub struct PyPoincare {
    pub poincare: Poincare,
}

#[pymethods]
impl PyPoincare {
    /// Creates a new [`PyPoincare`] object.
    #[new]
    pub fn new(init: PyPoincareInit, mapping: &PyMapping) -> Self {
        Self {
            poincare: Poincare::new(init.poincare_init, mapping.mapping),
        }
    }

    pub fn run(
        &mut self,
        qfactor: &PyQfactor,
        bfield: &PyBfield,
        current: &PyCurrent,
        per: &PyPerturbation,
    ) -> Result<(), PyPoincareError> {
        Ok(self.poincare.run(
            &qfactor.qfactor,
            &bfield.bfield,
            &current.current,
            &per.perturbation,
        )?)
    }

    // TODO: init and results getters

    #[getter]
    pub fn get_section(&self) -> String {
        format!("{:?}", self.poincare.mapping.section)
    }

    #[getter]
    pub fn get_alpha(&self) -> f64 {
        self.poincare.mapping.alpha
    }
}

to_numpy2D_impl!(PyPoincare, poincare, angles);
to_numpy2D_impl!(PyPoincare, poincare, fluxes);
