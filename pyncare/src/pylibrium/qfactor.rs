use std::path::PathBuf;

use equilibrium::Qfactor;
use rsl_interpolation::Accelerator;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::error::PyEqError;
use crate::{eval1D_impl, repr_impl, to_numpy1D_impl};

#[pyclass]
#[pyo3(name = "Qfactor")]
pub struct PyQfactor {
    pub qfactor: Qfactor,

    // Derived from [`Qfactor`].
    #[pyo3(get)]
    pub path: PathBuf,
    #[pyo3(get)]
    pub typ: String,
    #[pyo3(get)]
    pub psip_wall: f64,
    #[pyo3(get)]
    pub psi_wall: f64,

    // for Python-exposed evaluations
    pub acc: Accelerator,
}

#[pymethods]
impl PyQfactor {
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let qfactor = Qfactor::from_dataset(&path, typ)?;
        let typ = qfactor.typ.clone();
        let psip_wall = qfactor.psip_wall;
        let psi_wall = qfactor.psi_wall;

        Ok(Self {
            qfactor,
            path,
            typ,
            psip_wall,
            psi_wall,
            acc: Accelerator::new(),
        })
    }
}

repr_impl!(PyQfactor);
eval1D_impl!(PyQfactor, qfactor, q);
eval1D_impl!(PyQfactor, qfactor, psi);
to_numpy1D_impl!(PyQfactor, qfactor, q_data);
to_numpy1D_impl!(PyQfactor, qfactor, psip_data);
to_numpy1D_impl!(PyQfactor, qfactor, psi_data);
to_numpy1D_impl!(PyQfactor, qfactor, q_data_derived);

/// Remove duplicate entries
impl std::fmt::Debug for PyQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PyQfactor")
            .field("qfactor", &self.qfactor)
            .field("acc", &self.acc)
            .finish()
    }
}
