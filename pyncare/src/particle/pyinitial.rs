use particle::InitialConditions;
use utils::to_pyfloat_impl;

use pyo3::prelude::*;

#[derive(Debug)]
#[pyclass(name = "InitialConditions")]
pub struct PyInitialConditions {
    pub initial: InitialConditions,
    pub t0: f64,
    pub theta0: f64,
    pub psip0: f64,
    pub rho0: f64,
    pub zeta0: f64,
    pub mu: f64,
}

#[pymethods]
impl PyInitialConditions {
    #[new]
    #[rustfmt::skip]
    pub fn new(t0: f64, theta0: f64, psip0: f64, rho0: f64, zeta0: f64, mu: f64) -> Self {
        let initial = InitialConditions {t0, theta0, psip0, rho0, zeta0, mu};
        Self {initial, t0, theta0, psip0, rho0, zeta0, mu}
    }
}

to_pyfloat_impl!(PyInitialConditions, initial, t0);
to_pyfloat_impl!(PyInitialConditions, initial, theta0);
to_pyfloat_impl!(PyInitialConditions, initial, psip0);
to_pyfloat_impl!(PyInitialConditions, initial, rho0);
to_pyfloat_impl!(PyInitialConditions, initial, zeta0);
to_pyfloat_impl!(PyInitialConditions, initial, mu);
