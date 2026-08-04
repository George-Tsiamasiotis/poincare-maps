//! Defines the PyCurrent enum that holds one of the Current objects.

use numpy::{IntoPyArray, PyArray1};
use pyo3::{prelude::*, types::PyType};

use crate::*;

// ===============================================================================================

#[pyclass(name = "_PyLarCurrent", from_py_object, frozen, immutable_type)]
#[derive(Clone)]
pub struct PyLarCurrent(pub LarCurrent);

#[pyclass(name = "_PyNcCurrent", from_py_object, frozen, immutable_type)]
#[derive(Clone)]
pub struct PyNcCurrent(pub NcCurrent);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyCurrent", frozen, immutable_type, skip_from_py_object)]
pub enum PyCurrent {
    Lar(PyLarCurrent),
    Nc(PyNcCurrent),
}

#[pymethods] // Builders
impl PyCurrent {
    #[classmethod]
    pub fn build_lar<'py>(_: Bound<'py, PyType>) -> Result<Self> {
        Ok(Self::Lar(PyLarCurrent(LarCurrent::new())))
    }

    #[classmethod]
    pub fn build_nc<'py>(_: Bound<'py, PyType>, path: String, interp_type: String) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_1d_type(interp_type)?;
        let builder = NcCurrentBuilder::new(&path, typ);
        let current = builder.build()?;
        Ok(Self::Nc(PyNcCurrent(current)))
    }
}

/// References to the trait object and variants
impl PyCurrent {
    pub fn current(&self) -> &dyn Current {
        match self {
            PyCurrent::Lar(current) => &current.0,
            PyCurrent::Nc(current) => &current.0,
        }
    }

    pub fn lar(&self) -> Result<&LarCurrent> {
        match self {
            Self::Lar(current) => Ok(&current.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Current".into(),
                inner: "LarCurrent".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcCurrent> {
        match self {
            Self::Nc(current) => Ok(&current.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Current".into(),
                inner: "NcCurrent".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // EquilibriumObject Trait
impl PyCurrent {
    #[getter]
    pub fn object_type(&self) -> String {
        format!("{:?}", self.current().object_type())
    }

    #[getter]
    pub fn psi_state(&self) -> String {
        format!("{:?}", self.current().psi_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.current().psip_state())
    }
}

#[pymethods] // Current Trait
impl PyCurrent {
    pub fn g_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self.current().g_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn g_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self.current().g_of_psip(psip, &mut Accelerator::new())?)
    }

    pub fn i_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self.current().i_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn i_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self.current().i_of_psip(psip, &mut Accelerator::new())?)
    }

    pub fn dg_dpsi(&self, q: f64) -> Result<f64> {
        Ok(self.current().dg_dpsi(q, &mut Accelerator::new())?)
    }

    pub fn dg_dpsip(&self, q: f64) -> Result<f64> {
        Ok(self.current().dg_dpsip(q, &mut Accelerator::new())?)
    }

    pub fn di_dpsi(&self, psi: f64) -> Result<f64> {
        Ok(self.current().di_dpsi(psi, &mut Accelerator::new())?)
    }

    pub fn di_dpsip(&self, psip: f64) -> Result<f64> {
        Ok(self.current().di_dpsip(psip, &mut Accelerator::new())?)
    }
}

// ===============================================================================================

// #[pymethods] // Lar
// impl PyCurrent {}

#[pymethods] // Nc
impl PyCurrent {
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
    pub fn psi_last(&self) -> Result<f64> {
        self.nc()?.psi_last().ok_or(DexterError::AttributeError {
            obj: "NcCurrent".into(),
            attr: "psi_last".into(),
        })
    }

    #[getter]
    pub fn psip_last(&self) -> Result<f64> {
        self.nc()?.psip_last().ok_or(DexterError::AttributeError {
            obj: "NcCurrent".into(),
            attr: "psip_last".into(),
        })
    }

    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let current = self.nc()?;
        match name {
            "g_array" => return Ok(current.g_array().into_pyarray(py)),
            "i_array" => return Ok(current.i_array().into_pyarray(py)),
            "psi_array" => match current.psi_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            "psip_array" => match current.psip_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            _ => (),
        }
        Err(DexterError::AttributeError {
            obj: "NcCurrent".into(),
            attr: name.into(),
        })
    }
}

// ===============================================================================================

wrapper_debug_export!(PyLarCurrent);
wrapper_debug_export!(PyNcCurrent);

#[pymethods]
impl PyCurrent {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.current())
    }
}
