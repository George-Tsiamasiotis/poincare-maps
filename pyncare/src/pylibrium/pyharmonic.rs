use std::path::PathBuf;

use equilibrium::{Harmonic, HarmonicCache};
use rsl_interpolation::Accelerator;
use utils::{eval_harmonic_impl, repr_impl, to_numpy1D_impl};
use utils::{to_pyfloat_impl, to_pystr_impl};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::error::PyEqError;

#[derive(Clone)]
#[pyclass(name = "Harmonic")]
pub struct PyHarmonic {
    pub harmonic: Harmonic,
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

        Ok(Self {
            harmonic,
            cache: HarmonicCache::new(),
            acc: Accelerator::new(),
        })
    }
}

repr_impl!(PyHarmonic);
to_pystr_impl!(PyHarmonic, harmonic, typ);
to_pystr_impl!(PyHarmonic, harmonic, path);
to_pyfloat_impl!(PyHarmonic, harmonic, m);
to_pyfloat_impl!(PyHarmonic, harmonic, n);
to_pyfloat_impl!(PyHarmonic, harmonic, amax);
to_pyfloat_impl!(PyHarmonic, harmonic, phase);
to_pyfloat_impl!(PyHarmonic, harmonic, psip_wall);
to_numpy1D_impl!(PyHarmonic, harmonic, psip_data);
to_numpy1D_impl!(PyHarmonic, harmonic, a_data);
eval_harmonic_impl!(PyHarmonic, harmonic, h);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dpsip);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dtheta);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dzeta);
eval_harmonic_impl!(PyHarmonic, harmonic, dh_dt);

/// Remove Cache
impl std::fmt::Debug for PyHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PyHarmonic")
            .field("harmonic", &self.harmonic)
            .field("acc", &self.acc)
            .finish()
    }
}
