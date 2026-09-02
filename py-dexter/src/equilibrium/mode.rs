//! Defines the PyMode enum that holds one of the Mode objects.

use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use pyo3::{prelude::*, types::PyType};

use crate::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyFluteMode(Arc<FluteMode>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcFluteMode(Arc<NcFluteMode>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyMode", frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub enum PyMode {
    Flute(PyFluteMode),
    Nc(PyNcFluteMode),
}

#[pymethods] // Builders
impl PyMode {
    #[classmethod]
    pub fn build_flute<'py>(
        _: Bound<'py, PyType>,
        epsilon: f64,
        lcfs: &PyLastClosedFluxSurface,
        m: i64,
        n: i64,
        phase: f64,
    ) -> Result<Self> {
        let inner = PyFluteMode(Arc::new(FluteMode::new(epsilon, lcfs.0, m, n, phase)));
        Ok(Self::Flute(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(
        _: Bound<'py, PyType>,
        path: String,
        interp_type: String,
        m: i64,
        n: i64,
        phase_method: Bound<'py, PyAny>,
        analytical_threshold_index: usize,
    ) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_1d_type(interp_type)?;
        let phase_method = resolve_phase_method(phase_method)?;
        let builder = NcFluteModeBuilder::new(&path, typ, m, n)
            .with_phase_method(phase_method)
            .with_analytical_threshold_index(analytical_threshold_index);
        let mode = builder.build()?;
        let inner = PyNcFluteMode(Arc::new(mode));

        assert_eq!(Arc::strong_count(&inner.0), 1);
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyMode {
    pub fn mode(&self) -> &dyn Mode {
        match self {
            PyMode::Flute(mode) => mode.0.as_ref(),
            PyMode::Nc(mode) => mode.0.as_ref(),
        }
    }

    pub fn boxed_mode(&self) -> Box<dyn Mode> {
        match self {
            PyMode::Flute(mode) => Box::new(Arc::unwrap_or_clone(mode.0.clone())),
            PyMode::Nc(mode) => Box::new(Arc::unwrap_or_clone(mode.0.clone())),
        }
    }

    pub fn flute(&self) -> Result<&FluteMode> {
        match self {
            Self::Flute(mode) => Ok(&mode.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Mode".into(),
                inner: "FluteMode".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcFluteMode> {
        match self {
            Self::Nc(mode) => Ok(&mode.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Mode".into(),
                inner: "NcFluteMode".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // EquilibriumObject Trait
impl PyMode {
    #[getter]
    pub fn object_type(&self) -> String {
        format!("{:?}", self.mode().object_type())
    }

    #[getter]
    pub fn psi_state(&self) -> String {
        format!("{:?}", self.mode().psi_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.mode().psip_state())
    }
}

#[pymethods] // Mode Trait
#[rustfmt::skip]
impl PyMode {
    #[getter]
    pub fn psi_last(&self) -> Result<f64> {
        self.mode().psi_last().ok_or(DexterError::AttributeError {
            obj: "Mode".into(),
            attr: "psi_last".into(),
        })
    }

    #[getter]
    pub fn psip_last(&self) -> Result<f64> {
        self.mode().psip_last().ok_or(DexterError::AttributeError {
            obj: "Mode".into(),
            attr: "psip_last".into(),
        })
    }

    #[getter]
    pub fn m(&self) -> Result<i64> {
        Ok(self.mode().m())
    }

    #[getter]
    pub fn n(&self) -> Result<i64> {
        Ok(self.mode().n())
    }

    pub fn ampl_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().ampl_of_psi(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn ampl_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().ampl_of_psip(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn phase_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().phase_of_psi(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn phase_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().phase_of_psip(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn m_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().m_of_psi(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn m_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().m_of_psip(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_dpsi(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_dpsi(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_dpsip(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_dpsip(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_of_psi_dtheta(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_of_psi_dtheta(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_of_psip_dtheta(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_of_psip_dtheta(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_of_psi_dzeta(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_of_psi_dzeta(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_of_psip_dzeta(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_of_psip_dzeta(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_of_psi_dt(&self, psi: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_of_psi_dt(psi, theta, zeta, t, &mut self.mode().generate_cache())?)
    }

    pub fn dm_of_psip_dt(&self, psip: f64, theta: f64, zeta: f64, t: f64) -> Result<f64> {
        Ok(self.mode().dm_of_psip_dt(psip, theta, zeta, t, &mut self.mode().generate_cache())?)
    }
}

// ===============================================================================================

#[pymethods] // Flute
impl PyMode {
    #[getter]
    pub fn lcfs(&self) -> Result<PyLastClosedFluxSurface> {
        Ok(PyLastClosedFluxSurface(self.flute()?.lcfs()))
    }

    #[getter]
    pub fn epsilon(&self) -> Result<f64> {
        Ok(self.flute()?.epsilon())
    }

    #[getter]
    pub fn phase(&self) -> Result<f64> {
        Ok(self.flute()?.phase())
    }
}

#[pymethods] // NcFlute
impl PyMode {
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
    pub fn phase_method(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.phase_method()))
    }

    #[getter]
    pub fn analytical_threshold_index(&self) -> Result<usize> {
        Ok(self.nc()?.analytical_threshold_index())
    }

    #[getter]
    pub fn phase_average(&self) -> Result<f64> {
        self.nc()?
            .phase_average()
            .ok_or(DexterError::AttributeError {
                obj: "NcFluteMode".into(),
                attr: "phase_average".into(),
            })
    }

    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let mode = self.nc()?;
        match name {
            "alpha_array" => return Ok(mode.alpha_array().into_pyarray(py)),
            "phase_array" => return Ok(mode.phase_array().into_pyarray(py)),
            "psi_array" => match mode.psi_array() {
                Some(array) => return Ok(array.into_pyarray(py)),
                None => (),
            },
            "psip_array" => match mode.psip_array() {
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
}

// ===============================================================================================

wrapper_debug_export!(PyFluteMode);
wrapper_debug_export!(PyNcFluteMode);

#[pymethods]
impl PyMode {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.mode())
    }
}
