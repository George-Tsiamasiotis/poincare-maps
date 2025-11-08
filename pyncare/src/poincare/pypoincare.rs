use poincare::Poincare;
use safe_unwrap::safe_unwrap;
use utils::{py_debug_impl, py_get_numpy2D, py_repr_impl};

use crate::{PyBfield, PyCurrents, PyMappingParameters, PyParticle, PyPerturbation, PyQfactor};
use crate::{PyPoincareError, PyPoincareInit};

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

#[pyclass(name = "Poincare")]
pub struct PyPoincare(pub Poincare);

#[pymethods]
impl PyPoincare {
    /// Creates a new PyPoincare object.
    #[new]
    pub fn new(init: PyPoincareInit, params: &PyMappingParameters) -> Self {
        Self(Poincare::new(init.0, params.0))
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

    pub fn __getitem__(&self, n: usize) -> PyParticle {
        PyParticle(
            self.0
                .particles
                .get(n)
                .expect("Particle index out of bounds")
                .clone(),
        )
    }

    #[getter]
    pub fn get_section(&self) -> String {
        format!("{:?}", self.0.params.section)
    }

    #[getter]
    pub fn get_alpha(&self) -> f64 {
        self.0.params.alpha
    }

    #[getter]
    pub fn get_intersections(&self) -> usize {
        self.0.params.intersections
    }
}

py_debug_impl!(PyPoincare);
py_repr_impl!(PyPoincare);
py_get_numpy2D!(PyPoincare, angles);
py_get_numpy2D!(PyPoincare, fluxes);
