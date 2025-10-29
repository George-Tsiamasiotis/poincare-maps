use std::path::PathBuf;

use equilibrium::Bfield;
use rsl_interpolation::{Accelerator, Cache};

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;

use crate::error::PyEqError;
use crate::{eval2D_impl, repr_impl, to_numpy1D_impl, to_numpy2D_impl};

#[pyclass]
#[derive(Debug)]
#[pyo3(name = "Bfield")]
pub struct PyBfield {
    pub bfield: Bfield,

    // Derived from [`Bfield`].
    #[pyo3(get)]
    pub path: PathBuf,
    #[pyo3(get)]
    pub typ: String,
    #[pyo3(get)]
    pub baxis: f64,
    #[pyo3(get)]
    pub raxis: f64,
    #[pyo3(get)]
    pub psip_wall: f64,
    #[pyo3(get)]
    pub psi_wall: f64,

    // for Python-exposed evaluations
    pub xacc: Accelerator,
    pub yacc: Accelerator,
    pub cache: Cache<f64>,
}

#[pymethods]
impl PyBfield {
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let bfield = Bfield::from_dataset(&path, typ)?;
        let typ = bfield.typ.clone();
        let baxis = bfield.baxis;
        let raxis = bfield.baxis;
        let psip_wall = bfield.psip_wall;
        let psi_wall = bfield.psi_wall;

        Ok(Self {
            bfield,
            path,
            typ,
            baxis,
            raxis,
            psip_wall,
            psi_wall,
            xacc: Accelerator::new(),
            yacc: Accelerator::new(),
            cache: Cache::new(),
        })
    }
}

repr_impl!(PyBfield);
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
