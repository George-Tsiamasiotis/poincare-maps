use std::path::PathBuf;

use equilibrium::Bfield;
use rsl_interpolation::{Accelerator, Cache};
use utils::{eval2D_impl, repr_impl, to_numpy1D_impl, to_numpy2D_impl};
use utils::{to_pyfloat_impl, to_pystr_impl};

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;

use crate::error::PyEqError;

#[pyclass(name = "Bfield")]
pub struct PyBfield {
    pub bfield: Bfield,

    // for Python-exposed evaluations
    pub xacc: Accelerator,
    pub yacc: Accelerator,
    pub cache: Cache<f64>,
}

#[pymethods]
impl PyBfield {
    /// Creates a new PyBfield object.
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let bfield = Bfield::from_dataset(&path, typ)?;

        Ok(Self {
            bfield,
            xacc: Accelerator::new(),
            yacc: Accelerator::new(),
            cache: Cache::new(),
        })
    }
}

repr_impl!(PyBfield);
to_pystr_impl!(PyBfield, bfield, typ);
to_pystr_impl!(PyBfield, bfield, path);
to_pyfloat_impl!(PyBfield, bfield, baxis);
to_pyfloat_impl!(PyBfield, bfield, raxis);
to_pyfloat_impl!(PyBfield, bfield, psi_wall);
to_pyfloat_impl!(PyBfield, bfield, psip_wall);
eval2D_impl!(PyBfield, bfield, b);
eval2D_impl!(PyBfield, bfield, db_dtheta);
eval2D_impl!(PyBfield, bfield, db_dpsip);
eval2D_impl!(PyBfield, bfield, d2b_dtheta2);
eval2D_impl!(PyBfield, bfield, d2b_dpsip2);
eval2D_impl!(PyBfield, bfield, d2b_dpsip_dtheta);
to_numpy1D_impl!(PyBfield, bfield, psip_data);
to_numpy1D_impl!(PyBfield, bfield, theta_data);
to_numpy2D_impl!(PyBfield, bfield, b_data);
to_numpy2D_impl!(PyBfield, bfield, r_data);
to_numpy2D_impl!(PyBfield, bfield, z_data);
to_numpy2D_impl!(PyBfield, bfield, db_dpsip_data);
to_numpy2D_impl!(PyBfield, bfield, db_dtheta_data);

/// Remove Cache
impl std::fmt::Debug for PyBfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PyBfield")
            .field("bfield", &self.bfield)
            .field("xacc", &self.xacc)
            .field("yacc", &self.yacc)
            .finish()
    }
}
