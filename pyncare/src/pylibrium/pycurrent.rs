use std::path::PathBuf;

use equilibrium::Current;
use rsl_interpolation::Accelerator;
use utils::{eval1D_impl, repr_impl, to_numpy1D_impl};
use utils::{to_pyfloat_impl, to_pystr_impl};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::error::PyEqError;

#[derive(Debug)]
#[pyclass(name = "Current")]
pub struct PyCurrent {
    pub current: Current,
    // for Python-exposed evaluations
    pub acc: Accelerator,
}

#[pymethods]
impl PyCurrent {
    /// Creates a new PyCurrent object.
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let current = Current::from_dataset(&path, typ)?;

        Ok(Self {
            current,
            acc: Accelerator::new(),
        })
    }
}

repr_impl!(PyCurrent);
to_pystr_impl!(PyCurrent, current, typ);
to_pystr_impl!(PyCurrent, current, path);
to_pyfloat_impl!(PyCurrent, current, psip_wall);
eval1D_impl!(PyCurrent, current, g);
eval1D_impl!(PyCurrent, current, i);
eval1D_impl!(PyCurrent, current, dg_dpsip);
eval1D_impl!(PyCurrent, current, di_dpsip);
to_numpy1D_impl!(PyCurrent, current, psip_data);
to_numpy1D_impl!(PyCurrent, current, g_data);
to_numpy1D_impl!(PyCurrent, current, i_data);
