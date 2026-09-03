//! Defines the PyPerturbation container.

use std::sync::Arc;

use pyo3::prelude::*;
use pyo3::types::PyList;

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(
    name = "_PyPerturbation",
    frozen,
    immutable_type,
    sequence,
    from_py_object
)]
#[derive(Clone)]
pub struct PyPerturbation(Arc<Perturbation>);

#[pymethods]
impl PyPerturbation {
    #[new]
    pub fn new<'py>(modes: Bound<'py, PyList>) -> Result<Self> {
        let pymodes: Vec<PyMode> = modes.into_iter().map(|m| m.extract().unwrap()).collect();
        let modes: DynModes = pymodes.iter().map(|m| m.boxed_mode()).collect();
        let perturbation = Perturbation::new(modes);
        Ok(PyPerturbation(Arc::new(perturbation)))
    }

    pub fn __len__(&self) -> usize {
        self.0.count()
    }
}

// ===============================================================================================

#[pymethods] // Evaluations
impl PyPerturbation {
    pub fn p_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .p_of_psi(psi, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn p_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .p_of_psip(psip, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_dpsi(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_dpsi(psi, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_dpsip(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_dpsip(psip, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_of_psi_dtheta(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_of_psi_dtheta(psi, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_of_psip_dtheta(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_of_psip_dtheta(psip, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_of_psi_dzeta(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_of_psi_dzeta(psi, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_of_psip_dzeta(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_of_psip_dzeta(psip, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_of_psi_dt(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_of_psi_dt(psi, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn dp_of_psip_dt(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self
            .0
            .dp_of_psip_dt(psip, theta, zeta, t, &mut self.0.generate_caches())?)
    }
}

// ===============================================================================================

wrapper_debug_export!(PyPerturbation);

#[pymethods]
impl PyPerturbation {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.0)
    }
}
