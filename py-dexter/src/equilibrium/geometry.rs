//! Defines the PyGeometry enum that holds one of the Geometry objects.

use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::{prelude::*, types::PyType};

use crate::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyLarGeometry(Arc<LarGeometry>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcGeometry(Arc<NcGeometry>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyGeometry", frozen, immutable_type)]
pub enum PyGeometry {
    Lar(PyLarGeometry),
    Nc(PyNcGeometry),
}

#[pymethods] // Builders
impl PyGeometry {
    #[classmethod]
    pub fn build_lar<'py>(
        _: Bound<'py, PyType>,
        baxis: f64,
        raxis: f64,
        rlast: f64,
    ) -> Result<Self> {
        let inner = PyLarGeometry(Arc::new(LarGeometry::new(baxis, raxis, rlast)));
        Ok(Self::Lar(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(
        _: Bound<'py, PyType>,
        path: String,
        interp1d_type: String,
        interp2d_type: String,
    ) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ1 = resolve_interpolation_1d_type(interp1d_type)?;
        let typ2 = resolve_interpolation_2d_type(interp2d_type)?;
        let builder = NcGeometryBuilder::new(&path, typ1, typ2);
        let geometry = builder.build()?;
        let inner = PyNcGeometry(Arc::new(geometry));
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyGeometry {
    pub fn geometry(&self) -> &dyn Geometry {
        match self {
            PyGeometry::Lar(geometry) => geometry.0.as_ref(),
            PyGeometry::Nc(geometry) => geometry.0.as_ref(),
        }
    }

    pub fn fluxcommute(&self) -> Result<&dyn FluxCommute> {
        match self {
            PyGeometry::Lar(_) => Err(DexterError::EvalError(
                "LarGeometry does not support flux commutation".into(),
            )),
            PyGeometry::Nc(geometry) => Ok(geometry.0.as_ref()),
        }
    }

    pub fn lar(&self) -> Result<&LarGeometry> {
        match self {
            Self::Lar(geometry) => Ok(&geometry.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Geometry".into(),
                inner: "LarGeometry".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcGeometry> {
        match self {
            Self::Nc(geometry) => Ok(&geometry.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Geometry".into(),
                inner: "NcGeometry".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // EquilibriumObject Trait
impl PyGeometry {
    #[getter]
    pub fn object_type(&self) -> String {
        format!("{:?}", self.geometry().object_type())
    }

    #[getter]
    pub fn psi_state(&self) -> String {
        format!("{:?}", self.geometry().psi_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.geometry().psip_state())
    }
}

#[pymethods] // FluxCommute Trait
impl PyGeometry {
    pub fn psip_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self
            .fluxcommute()?
            .psip_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn psi_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self
            .fluxcommute()?
            .psi_of_psip(psip, &mut Accelerator::new())?)
    }
}

#[pymethods] // Geometry Trait
#[rustfmt::skip]
impl PyGeometry {
    #[getter]
    pub fn baxis(&self) -> f64 {
        self.geometry().baxis()
    }

    #[getter]
    pub fn raxis(&self) -> f64 {
        self.geometry().raxis()
    }

    #[getter]
    pub fn zaxis(&self) -> f64 {
        self.geometry().zaxis()
    }

    #[getter]
    pub fn rgeo(&self) -> f64 {
        self.geometry().rgeo()
    }

    #[getter]
    pub fn rlast(&self) -> f64 {
        self.geometry().rlast()
    }

    #[getter]
    pub fn psi_last(&self) -> Result<f64> {
        self.geometry().psi_last().ok_or(DexterError::AttributeError {
            obj: "Geometry".into(),
            attr: "psi_last".into(),
        })
    }

    #[getter]
    pub fn psip_last(&self) -> Result<f64> {
        self.geometry().psip_last().ok_or(DexterError::AttributeError {
            obj: "Geometry".into(),
            attr: "psip_last".into(),
        })
    }

    pub fn r_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self.geometry().r_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn r_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self.geometry().r_of_psip(psip, &mut Accelerator::new())?)
    }

    pub fn psi_of_r(&self, r: f64) -> Result<f64> {
        Ok(self.geometry().psi_of_r(r, &mut Accelerator::new())?)
    }

    pub fn psip_of_r(&self, r: f64) -> Result<f64> {
        Ok(self.geometry().psip_of_r(r, &mut Accelerator::new())?)
    }

    pub fn rlab_of_psi(&self, psi: f64, theta: f64) -> Result<f64> {
        Ok(self.geometry().rlab_of_psi(psi, theta, &mut Accelerator2d::new())?)
    }

    pub fn rlab_of_psip(&self, psip: f64, theta: f64) -> Result<f64> {
        Ok(self.geometry().rlab_of_psip(psip, theta, &mut Accelerator2d::new())?)
    }

    pub fn zlab_of_psi(&self, psi: f64, theta: f64) -> Result<f64> {
        Ok(self.geometry().zlab_of_psi(psi, theta, &mut Accelerator2d::new())?)
    }

    pub fn zlab_of_psip(&self, psip: f64, theta: f64) -> Result<f64> {
        Ok(self.geometry().zlab_of_psip(psip, theta, &mut Accelerator2d::new())?)
    }

    pub fn jacobian_of_psi(&self, psi: f64, theta: f64) -> Result<f64> {
        Ok(self.geometry().jacobian_of_psi(psi, theta, &mut Accelerator2d::new())?)
    }

    pub fn jacobian_of_psip(&self, psip: f64, theta: f64) -> Result<f64> {
        Ok(self.geometry().jacobian_of_psip(psip, theta, &mut Accelerator2d::new())?)
    }

    #[getter]
    pub fn rlab_last<'py>(&self, py: Python<'py>)  -> Result<Bound<'py, PyArray1<f64>>> {
        Ok(self.geometry().rlab_last().into_pyarray(py))
    }

    #[getter]
    pub fn zlab_last<'py>(&self, py: Python<'py>)  -> Result<Bound<'py, PyArray1<f64>>> {
        Ok(self.geometry().zlab_last().into_pyarray(py))
    }
}

// ===============================================================================================

// #[pymethods] // Lar
// impl PyGeometry {}

#[pymethods] // Nc
impl PyGeometry {
    #[getter]
    pub fn path(&self) -> Result<String> {
        Ok(self.nc()?.path().to_str().unwrap_or_default().to_string())
    }

    #[getter]
    pub fn netcdf_version(&self) -> Result<String> {
        Ok(self.nc()?.netcdf_version().to_string())
    }

    #[getter]
    pub fn interp1d_type(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.interp1d_type()))
    }

    #[getter]
    pub fn interp2d_type(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.interp2d_type()))
    }

    #[getter]
    pub fn shape(&self) -> Result<(usize, usize)> {
        Ok(self.nc()?.shape())
    }

    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let geometry = self.nc()?;
        match name {
            "theta_array" => return Ok(geometry.theta_array().into_pyarray(py)),
            "r_array" => return Ok(geometry.r_array().into_pyarray(py)),
            "psi_array" => match geometry.psi_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            "psip_array" => match geometry.psip_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            _ => (),
        }
        Err(DexterError::AttributeError {
            obj: "NcGeometry".into(),
            attr: name.into(),
        })
    }

    pub fn get_array2d<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Bound<'py, PyArray2<f64>>> {
        let geometry = self.nc()?;
        match name {
            "rlab_array" => return Ok(geometry.rlab_array().into_pyarray(py)),
            "zlab_array" => return Ok(geometry.zlab_array().into_pyarray(py)),
            "jacobian_array" => return Ok(geometry.jacobian_array().into_pyarray(py)),
            _ => (),
        }
        Err(DexterError::AttributeError {
            obj: "NcGeometry".into(),
            attr: name.into(),
        })
    }
}

// ===============================================================================================

wrapper_debug_export!(PyLarGeometry);
wrapper_debug_export!(PyNcGeometry);

#[pymethods]
impl PyGeometry {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.geometry())
    }
}
