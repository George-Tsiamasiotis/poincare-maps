//! Defines the PyBfield enum that holds one of the Bfield objects.

use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::{prelude::*, types::PyType};

use crate::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyLarBfield(Arc<LarBfield>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcBfield(Arc<NcBfield>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyBfield", frozen, immutable_type)]
pub enum PyBfield {
    Lar(PyLarBfield),
    Nc(PyNcBfield),
}

#[pymethods] // Builders
impl PyBfield {
    #[classmethod]
    pub fn build_lar<'py>(_: Bound<'py, PyType>) -> Result<Self> {
        let inner = PyLarBfield(Arc::new(LarBfield::new()));
        Ok(Self::Lar(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(
        _: Bound<'py, PyType>,
        path: String,
        interp_type: String,
        padding: usize,
    ) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_2d_type(interp_type)?;
        let builder = NcBfieldBuilder::new(&path, typ).with_padding(padding);
        let bfield = builder.build()?;
        let inner = PyNcBfield(Arc::new(bfield));
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyBfield {
    pub fn bfield(&self) -> &dyn Bfield {
        match self {
            PyBfield::Lar(bfield) => bfield.0.as_ref(),
            PyBfield::Nc(bfield) => bfield.0.as_ref(),
        }
    }

    pub fn lar(&self) -> Result<&LarBfield> {
        match self {
            Self::Lar(bfield) => Ok(&bfield.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Bfield".into(),
                inner: "LarBfield".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcBfield> {
        match self {
            Self::Nc(bfield) => Ok(&bfield.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Bfield".into(),
                inner: "NcBfield".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // EquilibriumObject Trait
impl PyBfield {
    #[getter]
    pub fn object_type(&self) -> String {
        format!("{:?}", self.bfield().object_type())
    }

    #[getter]
    pub fn psi_state(&self) -> String {
        format!("{:?}", self.bfield().psi_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.bfield().psip_state())
    }
}

#[pymethods] // Bfield Trait
#[rustfmt::skip]
impl PyBfield {
    pub fn b_of_psi(&self, psi: f64, theta: f64) -> Result<f64> {
        Ok(self.bfield().b_of_psi(psi, theta, &mut Accelerator2d::new())?)
    }

    pub fn b_of_psip(&self, psip: f64, theta: f64) -> Result<f64> {
        Ok(self.bfield().b_of_psip(psip, theta, &mut Accelerator2d::new())?)
    }

    pub fn db_dpsi(&self, psi: f64, theta: f64) -> Result<f64> {
        Ok(self.bfield().db_dpsi(psi, theta, &mut Accelerator2d::new())?)
    }

    pub fn db_dpsip(&self, psip: f64, theta: f64) -> Result<f64> {
        Ok(self.bfield().db_dpsip(psip, theta, &mut Accelerator2d::new())?)
    }

    pub fn db_of_psi_dtheta(&self, psi: f64, theta: f64) -> Result<f64> {
        Ok(self.bfield().db_of_psi_dtheta(psi, theta, &mut Accelerator2d::new())?)
    }

    pub fn db_of_psip_dtheta(&self, psip: f64, theta: f64) -> Result<f64> {
        Ok(self.bfield().db_of_psip_dtheta(psip, theta, &mut Accelerator2d::new())?)
    }
}

// ===============================================================================================

// #[pymethods] // Lar
// impl PyBfield {}

#[pymethods] // Nc
impl PyBfield {
    #[getter]
    pub fn path(&self) -> Result<String> {
        Ok(self.nc()?.path().to_str().unwrap_or_default().to_string())
    }

    #[getter]
    pub fn netcdf_version(&self) -> Result<String> {
        Ok(self.nc()?.netcdf_version().to_string())
    }

    #[getter]
    pub fn interp_type(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.interp_type()))
    }

    #[getter]
    pub fn baxis(&self) -> Result<f64> {
        Ok(self.nc()?.baxis())
    }

    #[getter]
    pub fn padding(&self) -> Result<usize> {
        Ok(self.nc()?.padding())
    }

    #[getter]
    pub fn shape(&self) -> Result<(usize, usize)> {
        Ok(self.nc()?.shape())
    }

    #[getter]
    pub fn shape_padded(&self) -> Result<(usize, usize)> {
        Ok(self.nc()?.shape_padded())
    }

    #[getter]
    pub fn psi_last(&self) -> Result<f64> {
        self.nc()?.psi_last().ok_or(DexterError::AttributeError {
            obj: "NcBfield".into(),
            attr: "psi_last".into(),
        })
    }

    #[getter]
    pub fn psip_last(&self) -> Result<f64> {
        self.nc()?.psip_last().ok_or(DexterError::AttributeError {
            obj: "NcBfield".into(),
            attr: "psip_last".into(),
        })
    }

    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let bfield = self.nc()?;
        match name {
            "theta_array" => return Ok(bfield.theta_array().into_pyarray(py)),
            "theta_array_padded" => return Ok(bfield.theta_array_padded().into_pyarray(py)),
            "psi_array" => match bfield.psi_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            "psip_array" => match bfield.psip_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            _ => (),
        }
        Err(DexterError::AttributeError {
            obj: "NcBfield".into(),
            attr: name.into(),
        })
    }

    pub fn get_array2d<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Bound<'py, PyArray2<f64>>> {
        let bfield = self.nc()?;
        match name {
            "b_array" => return Ok(bfield.b_array().into_pyarray(py)),
            "b_array_padded" => return Ok(bfield.b_array_padded().into_pyarray(py)),
            _ => (),
        }
        Err(DexterError::AttributeError {
            obj: "NcBfield".into(),
            attr: name.into(),
        })
    }
}
