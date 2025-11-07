use poincare::Poincare;
use safe_unwrap::safe_unwrap;
use utils::py_get_numpy2D;

use crate::{PyBfield, PyCurrents, PyMappingParameters, PyPerturbation, PyQfactor};
use crate::{PyPoincareError, PyPoincareInit};

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

#[pyclass(name = "Poincare")]
pub struct PyPoincare(pub Poincare);

#[pymethods]
impl PyPoincare {
    /// Creates a new PyPoincare object.
    #[new]
    pub fn new(init: PyPoincareInit, mapping: &PyMappingParameters) -> Self {
        Self(Poincare::new(init.0, mapping.0))
    }

    pub fn run(
        &mut self,
        qfactor: &PyQfactor,
        bfield: &PyBfield,
        currents: &PyCurrents,
        perturbation: &PyPerturbation,
    ) -> Result<(), PyPoincareError> {
        Ok(self
            .0
            .run(&qfactor.0, &bfield.0, &currents.0, &perturbation.0)?)
    }

    #[getter]
    pub fn get_init(&self) -> PyPoincareInit {
        let init = &self.0.init;
        safe_unwrap!(
            "cannot fail if self.0.init exists",
            PyPoincareInit::new(
                init.thetas().to_vec(),
                init.psips().to_vec(),
                init.rhos().to_vec(),
                init.zetas().to_vec(),
                init.mus().to_vec(),
            )
        )
    }

    #[getter]
    pub fn get_section(&self) -> String {
        format!("{:?}", self.0.mapping.section)
    }

    #[getter]
    pub fn get_alpha(&self) -> f64 {
        self.0.mapping.alpha
    }

    #[getter]
    pub fn get_intersections(&self) -> usize {
        self.0.mapping.intersections
    }
}

// TODO:
// py_debug_impl!(PyPoincare);
// py_repr_impl!(PyPoincare);
py_get_numpy2D!(PyPoincare, angles);
py_get_numpy2D!(PyPoincare, fluxes);
