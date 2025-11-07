use equilibrium::Qfactor;
use rsl_interpolation::Accelerator;
use safe_unwrap::safe_unwrap;
use utils::{py_eval1D, py_get_numpy1D, py_get_numpy1D_fallible, py_repr_impl};
use utils::{py_get_float, py_get_path, py_get_typ};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::PyEqError;

#[pyclass(name = "Qfactor")]
pub struct PyQfactor(pub Qfactor);

#[pymethods]
impl PyQfactor {
    /// Creates a new PyQFactor wrapper object.
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Qfactor::from_dataset(&path, typ)?))
    }
}

py_repr_impl!(PyQfactor);
py_get_path!(PyQfactor);
py_get_typ!(PyQfactor);
py_get_float!(PyQfactor, psip_wall);
py_get_float!(PyQfactor, psi_wall);
py_eval1D!(PyQfactor, q);
py_eval1D!(PyQfactor, psi);
py_get_numpy1D!(PyQfactor, psip_data);
py_get_numpy1D!(PyQfactor, psi_data);
py_get_numpy1D!(PyQfactor, q_data);
py_get_numpy1D_fallible!(PyQfactor, q_data_derived);
