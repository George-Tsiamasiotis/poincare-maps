use poincare::PoincareInit;
use utils::{repr_impl, to_numpy1D_impl};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::PyPoincareError;

#[pyclass(name = "PoincareInit")]
pub struct PyPoincareInit {
    poincare_init: PoincareInit,
}

#[pymethods]
impl PyPoincareInit {
    /// Creates a new PyPoincareInit object.
    #[new]
    pub fn new(
        thetas: Vec<f64>,
        psips: Vec<f64>,
        rhos: Vec<f64>,
        zetas: Vec<f64>,
        mus: Vec<f64>,
    ) -> Result<Self, PyPoincareError> {
        Ok(Self {
            poincare_init: PoincareInit::build(&thetas, &psips, &rhos, &zetas, &mus)?,
        })
    }
}

impl std::fmt::Debug for PyPoincareInit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.poincare_init.fmt(f)
    }
}

repr_impl!(PyPoincareInit);
to_numpy1D_impl!(PyPoincareInit, poincare_init, thetas);
to_numpy1D_impl!(PyPoincareInit, poincare_init, psips);
to_numpy1D_impl!(PyPoincareInit, poincare_init, rhos);
to_numpy1D_impl!(PyPoincareInit, poincare_init, zetas);
to_numpy1D_impl!(PyPoincareInit, poincare_init, mus);
