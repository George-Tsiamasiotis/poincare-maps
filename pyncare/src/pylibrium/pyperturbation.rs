use equilibrium::{Harmonic, HarmonicCache, Perturbation};
use rsl_interpolation::Accelerator;
use utils::{py_eval_perturbation, py_repr_impl};

use pyo3::prelude::*;
use pyo3::types::PyList;

use crate::PyEqError;
use crate::PyHarmonic;

#[pyclass(name = "Perturbation")]
pub struct PyPerturbation(Perturbation);

#[pymethods]
impl PyPerturbation {
    /// Creates a new PyPerturbation wrapper object.
    #[new]
    pub fn new_py<'py>(pyharmonics: Bound<'py, PyList>) -> Result<Self, PyEqError> {
        let pyharmonics_vec: Vec<PyHarmonic> = pyharmonics
            .iter()
            .map(|ph| {
                ph.extract()
                    .expect("Could not extract 'PyHarmonic' from python list")
            })
            .collect();
        let harmonics_vec: Vec<Harmonic> = pyharmonics_vec
            .clone()
            .iter()
            .map(|ph| ph.0.clone())
            .collect();

        Ok(Self(Perturbation::from_harmonics(&harmonics_vec)))
    }

    /// Makes PyPerturbation indexable
    pub fn __getitem__(&self, index: usize) -> PyHarmonic {
        PyHarmonic::from(
            self.0
                .harmonics
                .get(index)
                .expect("Harmonic index out of bounds"),
        )
    }

    /// Returns the number of Harmonics.
    pub fn len(&self) -> usize {
        self.0.harmonics.len()
    }

    /// Return `true` if the Perturbation contains no Harmonics.
    pub fn is_empty(&self) -> bool {
        self.0.harmonics.is_empty()
    }
}

py_repr_impl!(PyPerturbation);
py_eval_perturbation!(PyPerturbation, p);
py_eval_perturbation!(PyPerturbation, dp_dpsip);
py_eval_perturbation!(PyPerturbation, dp_dtheta);
py_eval_perturbation!(PyPerturbation, dp_dzeta);
py_eval_perturbation!(PyPerturbation, dp_dt);
