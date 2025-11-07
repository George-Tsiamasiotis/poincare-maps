use equilibrium::Currents;
use rsl_interpolation::Accelerator;
use safe_unwrap::safe_unwrap;
use utils::{py_eval1D, py_get_numpy1D, py_repr_impl};
use utils::{py_get_float, py_get_path, py_get_typ};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::PyEqError;

#[pyclass(name = "Current")]
pub struct PyCurrents(pub Currents);

#[pymethods]
impl PyCurrents {
    /// Creates a new PyCurrents wrapper object.
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Currents::from_dataset(&path, typ)?))
    }
}

py_repr_impl!(PyCurrents);
py_get_typ!(PyCurrents);
py_get_path!(PyCurrents);
py_get_float!(PyCurrents, psip_wall);
py_eval1D!(PyCurrents, g);
py_eval1D!(PyCurrents, i);
py_eval1D!(PyCurrents, dg_dpsip);
py_eval1D!(PyCurrents, di_dpsip);
py_get_numpy1D!(PyCurrents, psip_data);
py_get_numpy1D!(PyCurrents, g_data);
py_get_numpy1D!(PyCurrents, i_data);
