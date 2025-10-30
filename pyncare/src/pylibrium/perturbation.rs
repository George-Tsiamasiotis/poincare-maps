use equilibrium::{Harmonic, HarmonicCache, Perturbation};
use pyo3::types::PyList;
use rsl_interpolation::Accelerator;

use crate::error::PyEqError;
use crate::{eval_harmonic_impl, repr_impl, PyHarmonic};

use pyo3::prelude::*;

#[pyclass]
#[pyo3(name = "Perturbation")]
pub struct PyPerturbation {
    pub perturbation: Perturbation,
    #[pyo3(get)]
    pub harmonics: Vec<PyHarmonic>,

    // for Python-exposed evaluations
    pub cache: Vec<HarmonicCache>,
    pub acc: Accelerator,
}

#[pymethods]
impl PyPerturbation {
    #[new]
    pub fn new_py<'py>(harmonics: Bound<'py, PyList>) -> Result<Self, PyEqError> {
        let pyharmonics_vec: Vec<PyHarmonic> =
            harmonics.iter().map(|h| h.extract().unwrap()).collect();
        let harmonics_vec: Vec<Harmonic> = pyharmonics_vec
            .clone()
            .iter()
            .map(|p| p.harmonic.clone())
            .collect();

        Ok(Self {
            perturbation: Perturbation::from_harmonics(&harmonics_vec),
            harmonics: pyharmonics_vec,
            acc: Accelerator::new(),
            cache: vec![HarmonicCache::new(); harmonics_vec.len()],
        })
    }

    /// Makes PyPerturbation indexable
    pub fn __getitem__(&self, key: usize) -> PyHarmonic {
        self.harmonics
            .get(key)
            .expect("Out of bounds access")
            .clone()
    }
}

repr_impl!(PyPerturbation);
eval_harmonic_impl!(PyPerturbation, perturbation, p);
eval_harmonic_impl!(PyPerturbation, perturbation, dp_dpsip);
eval_harmonic_impl!(PyPerturbation, perturbation, dp_dtheta);
eval_harmonic_impl!(PyPerturbation, perturbation, dp_dzeta);
eval_harmonic_impl!(PyPerturbation, perturbation, dp_dt);

impl std::fmt::Debug for PyPerturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.harmonics.fmt(f)
    }
}
