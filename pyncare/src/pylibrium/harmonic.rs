use equilibrium::{Harmonic, HarmonicCache};
use rsl_interpolation::Accelerator;
use std::path::PathBuf;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::error::PyEqError;
use crate::{eval_harmonic_impl, repr_impl, to_numpy1D_impl};

#[pyclass]
#[pyo3(name = "Harmonic")]
pub struct PyHarmonic {
    pub harmonic: Harmonic,

    // Derived from [`Harmonic`].
    #[pyo3(get)]
    pub path: PathBuf,
    #[pyo3(get)]
    pub typ: String,
    #[pyo3(get)]
    pub m: f64,
    #[pyo3(get)]
    pub n: f64,
    #[pyo3(get)]
    pub phase: f64,
    #[pyo3(get)]
    pub amax: f64,
    #[pyo3(get)]
    pub psip_wall: f64,

    // for Python-exposed evaluations
    pub cache: HarmonicCache,
    pub acc: Accelerator,
}

#[pymethods]
impl PyHarmonic {
    #[new]
    pub fn new(path: &str, typ: &str, m: f64, n: f64, phase: f64) -> Result<Self, PyEqError> {
        let path = PathBuf::from(path);
        let harmonic = Harmonic::from_dataset(&path, typ, m, n, phase)?;
        let amax = harmonic.amax;
        let psip_wall = harmonic.psip_wall;

        Ok(Self {
            harmonic,
            path,
            typ: typ.into(),
            m,
            n,
            phase,
            amax,
            psip_wall,
            cache: HarmonicCache::new(),
            acc: Accelerator::new(),
        })
    }
}

repr_impl!(PyHarmonic);
to_numpy1D_impl!(PyHarmonic, harmonic, psip_data);
to_numpy1D_impl!(PyHarmonic, harmonic, a_data);
eval_harmonic_impl!(PyHarmonic, harmonic, h);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dpsip);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dtheta);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dzeta);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dt);

/// Remove duplicate entries
impl std::fmt::Debug for PyHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PyHarmonic")
            .field("harmonic", &self.harmonic)
            .field("acc", &self.acc)
            .finish()
    }
}
