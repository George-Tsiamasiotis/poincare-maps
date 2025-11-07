use equilibrium::Bfield;
use rsl_interpolation::{Accelerator, Cache};
use safe_unwrap::safe_unwrap;
use utils::{py_eval2D, py_get_numpy1D, py_get_numpy2D, py_get_numpy2D_fallible, py_repr_impl};
use utils::{py_get_float, py_get_path, py_get_typ};

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;

use crate::PyEqError;

#[pyclass(name = "Bfield")]
pub struct PyBfield(pub Bfield);

#[pymethods]
impl PyBfield {
    /// Creates a new PyBfield wrapper object.
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Bfield::from_dataset(&path, typ)?))
    }
}

py_repr_impl!(PyBfield);
py_get_typ!(PyBfield);
py_get_path!(PyBfield);
py_get_float!(PyBfield, psip_wall);
py_get_float!(PyBfield, baxis);
py_get_float!(PyBfield, raxis);
py_eval2D!(PyBfield, b);
py_eval2D!(PyBfield, db_dtheta);
py_eval2D!(PyBfield, db_dpsip);
py_eval2D!(PyBfield, d2b_dtheta2);
py_eval2D!(PyBfield, d2b_dpsip2);
py_eval2D!(PyBfield, d2b_dpsip_dtheta);
py_get_numpy1D!(PyBfield, psip_data);
py_get_numpy1D!(PyBfield, theta_data);
py_get_numpy2D!(PyBfield, b_data);
py_get_numpy2D!(PyBfield, r_data);
py_get_numpy2D!(PyBfield, z_data);
py_get_numpy2D_fallible!(PyBfield, db_dpsip_data);
py_get_numpy2D_fallible!(PyBfield, db_dtheta_data);
