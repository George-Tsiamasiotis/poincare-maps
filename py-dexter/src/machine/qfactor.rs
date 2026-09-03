//! Defines the PyQfactor enum that holds one of the Qfactor objects.

use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use pyo3::{prelude::*, types::PyType};

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyUnityQfactor(Arc<UnityQfactor>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyParabolicQfactor(Arc<ParabolicQfactor>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcQfactor(Arc<NcQfactor>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyQfactor", frozen, immutable_type)]
pub enum PyQfactor {
    Unity(PyUnityQfactor),
    Parabolic(PyParabolicQfactor),
    Nc(PyNcQfactor),
}

#[pymethods] // Builders
impl PyQfactor {
    #[classmethod]
    pub fn build_unity<'py>(_: Bound<'py, PyType>, lcfs: &PyLastClosedFluxSurface) -> Result<Self> {
        let inner = PyUnityQfactor(Arc::new(UnityQfactor::new(lcfs.0)));
        Ok(Self::Unity(inner))
    }

    #[classmethod]
    pub fn build_parabolic<'py>(
        _: Bound<'py, PyType>,
        qaxis: f64,
        qlast: f64,
        lcfs: &PyLastClosedFluxSurface,
    ) -> Result<Self> {
        let inner = PyParabolicQfactor(Arc::new(ParabolicQfactor::new(qaxis, qlast, lcfs.0)));
        Ok(Self::Parabolic(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(_: Bound<'py, PyType>, path: String, interp_type: String) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_1d_type(interp_type)?;
        let builder = NcQfactorBuilder::new(&path, typ);
        let qfactor = builder.build()?;
        let inner = PyNcQfactor(Arc::new(qfactor));
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyQfactor {
    pub fn inner(&self) -> &dyn Qfactor {
        match self {
            PyQfactor::Unity(qfactor) => qfactor.0.as_ref(),
            PyQfactor::Parabolic(qfactor) => qfactor.0.as_ref(),
            PyQfactor::Nc(qfactor) => qfactor.0.as_ref(),
        }
    }

    pub fn unity(&self) -> Result<&UnityQfactor> {
        match self {
            Self::Unity(qfactor) => Ok(&qfactor.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Qfactor".into(),
                inner: "UnityQfactor".into(),
            }),
        }
    }

    pub fn parabolic(&self) -> Result<&ParabolicQfactor> {
        match self {
            Self::Parabolic(qfactor) => Ok(&qfactor.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Qfactor".into(),
                inner: "ParabolicQfactor".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcQfactor> {
        match self {
            Self::Nc(qfactor) => Ok(&qfactor.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Qfactor".into(),
                inner: "NcQfactor".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // EquilibriumObject Trait
impl PyQfactor {
    #[getter]
    pub fn machine_type(&self) -> String {
        format!("{:?}", self.inner().machine_type())
    }

    #[getter]
    pub fn psi_state(&self) -> String {
        format!("{:?}", self.inner().psi_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.inner().psip_state())
    }
}

#[pymethods] // FluxCommute Trait
impl PyQfactor {
    pub fn psip_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self.inner().psip_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn psi_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self.inner().psi_of_psip(psip, &mut Accelerator::new())?)
    }
}

#[pymethods] // Qfactor Trait
impl PyQfactor {
    #[getter]
    pub fn psi_last(&self) -> f64 {
        self.inner().psi_last()
    }

    #[getter]
    pub fn psip_last(&self) -> f64 {
        self.inner().psip_last()
    }

    #[getter]
    pub fn qlast(&self) -> f64 {
        self.inner().qlast()
    }

    #[getter]
    pub fn qaxis(&self) -> f64 {
        self.inner().qaxis()
    }

    pub fn q_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self.inner().q_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn q_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self.inner().q_of_psip(psip, &mut Accelerator::new())?)
    }

    pub fn dpsip_dpsi(&self, psi: f64) -> Result<f64> {
        Ok(self.inner().dpsip_dpsi(psi, &mut Accelerator::new())?)
    }

    pub fn dpsi_dpsip(&self, psip: f64) -> Result<f64> {
        Ok(self.inner().dpsi_dpsip(psip, &mut Accelerator::new())?)
    }

    pub fn psi_of_q(&self, q: f64) -> Result<f64> {
        Ok(self.inner().psi_of_q(q, &mut Accelerator::new())?)
    }

    pub fn psip_of_q(&self, q: f64) -> Result<f64> {
        Ok(self.inner().psip_of_q(q, &mut Accelerator::new())?)
    }

    pub fn iota_of_psi(&self, psi: f64) -> Result<f64> {
        Ok(self.inner().iota_of_psi(psi, &mut Accelerator::new())?)
    }

    pub fn iota_of_psip(&self, psip: f64) -> Result<f64> {
        Ok(self.inner().iota_of_psip(psip, &mut Accelerator::new())?)
    }
}

// ===============================================================================================

// #[pymethods] // Unity
// impl PyQfactor {}
//
// #[pymethods] // Parabolic
// impl PyQfactor {}

#[pymethods] // Nc
impl PyQfactor {
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

    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let qfactor = self.nc()?;
        match name {
            "q_array" => return Ok(qfactor.q_array().into_pyarray(py)),
            "psi_array" => match qfactor.psi_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            "psip_array" => match qfactor.psip_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            _ => (),
        }
        Err(DexterError::AttributeError {
            obj: "NcQfactor".into(),
            attr: name.into(),
        })
    }
}

// ===============================================================================================

wrapper_debug_export!(PyUnityQfactor);
wrapper_debug_export!(PyParabolicQfactor);
wrapper_debug_export!(PyNcQfactor);

#[pymethods]
impl PyQfactor {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.inner())
    }
}
