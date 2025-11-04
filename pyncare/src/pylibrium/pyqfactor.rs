use std::path::PathBuf;

use equilibrium::Qfactor;
use rsl_interpolation::Accelerator;
use utils::{eval1D_impl, repr_impl, to_numpy1D_impl};
use utils::{to_pyfloat_impl, to_pystr_impl};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::error::PyEqError;

#[derive(Debug)]
#[pyclass(name = "Qfactor")]
pub struct PyQfactor {
    pub qfactor: Qfactor,
    // for Python-exposed evaluations
    pub acc: Accelerator,
}

#[pymethods]
impl PyQfactor {
    /// Creates a new PyQFactor object.
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let qfactor = Qfactor::from_dataset(&path, typ)?;

        Ok(Self {
            qfactor,
            acc: Accelerator::new(),
        })
    }
}

repr_impl!(PyQfactor);
to_pystr_impl!(PyQfactor, qfactor, typ);
to_pystr_impl!(PyQfactor, qfactor, path);
to_pyfloat_impl!(PyQfactor, qfactor, psi_wall);
to_pyfloat_impl!(PyQfactor, qfactor, psip_wall);
eval1D_impl!(PyQfactor, qfactor, q);
eval1D_impl!(PyQfactor, qfactor, psi);
to_numpy1D_impl!(PyQfactor, qfactor, q_data);
to_numpy1D_impl!(PyQfactor, qfactor, psip_data);
to_numpy1D_impl!(PyQfactor, qfactor, psi_data);
to_numpy1D_impl!(PyQfactor, qfactor, q_data_derived);
