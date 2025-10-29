use std::path::PathBuf;

use equilibrium::Current;
use rsl_interpolation::Accelerator;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::error::PyEqError;
use crate::{eval1D_impl, repr_impl, to_numpy1D_impl};

#[pyclass]
#[derive(Debug)]
#[pyo3(name = "Current")]
pub struct PyCurrent {
    pub current: Current,

    // Derived from [`Current`].
    #[pyo3(get)]
    pub path: PathBuf,
    #[pyo3(get)]
    pub typ: String,
    #[pyo3(get)]
    pub psip_wall: f64,

    // for Python-exposed evaluations
    pub acc: Accelerator,
}

#[pymethods]
impl PyCurrent {
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let current = Current::from_dataset(&path, typ)?;
        let typ = current.typ.clone();
        let psip_wall = current.psip_wall;

        Ok(Self {
            current,
            path,
            typ,
            psip_wall,
            acc: Accelerator::new(),
        })
    }
}

repr_impl!(PyCurrent);
eval1D_impl!(PyCurrent, current, g);
eval1D_impl!(PyCurrent, current, i);
eval1D_impl!(PyCurrent, current, dg_dpsip);
eval1D_impl!(PyCurrent, current, di_dpsip);
to_numpy1D_impl!(PyCurrent, current, psip_data);
to_numpy1D_impl!(PyCurrent, current, g_data);
to_numpy1D_impl!(PyCurrent, current, i_data);
to_numpy1D_impl!(PyCurrent, current, dg_dpsip_data);
to_numpy1D_impl!(PyCurrent, current, di_dpsip_data);
