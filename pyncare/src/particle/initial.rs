use particle::InitialConditions;

use pyo3::prelude::*;

#[pyclass]
#[derive(Debug)]
#[pyo3(name = "InitialConditions")]
pub struct PyInitialConditions {
    pub initial: InitialConditions,
    #[pyo3(get)]
    pub t0: f64,
    #[pyo3(get)]
    pub theta0: f64,
    #[pyo3(get)]
    pub psip0: f64,
    #[pyo3(get)]
    pub rho0: f64,
    #[pyo3(get)]
    pub zeta0: f64,
    #[pyo3(get)]
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
