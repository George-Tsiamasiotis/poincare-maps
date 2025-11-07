use poincare::PoincareInit;
use utils::{py_debug_impl, py_get_numpy1D, py_repr_impl};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::PyPoincareError;

#[derive(Clone)]
#[pyclass(name = "PoincareInit")]
pub struct PyPoincareInit(pub PoincareInit);

#[pymethods]
impl PyPoincareInit {
    /// Creates a new PyPoincareInit object.
    ///
    /// Use Vec<f64> and let pyo3 do the rest.
    #[new]
    pub fn new(
        thetas: Vec<f64>,
        psips: Vec<f64>,
        rhos: Vec<f64>,
        zetas: Vec<f64>,
        mus: Vec<f64>,
    ) -> Result<Self, PyPoincareError> {
        Ok(Self(PoincareInit::build(
            &thetas, &psips, &rhos, &zetas, &mus,
        )?))
    }
}

py_debug_impl!(PyPoincareInit);
py_repr_impl!(PyPoincareInit);
py_get_numpy1D!(PyPoincareInit, thetas);
py_get_numpy1D!(PyPoincareInit, psips);
py_get_numpy1D!(PyPoincareInit, rhos);
py_get_numpy1D!(PyPoincareInit, zetas);
py_get_numpy1D!(PyPoincareInit, mus);
