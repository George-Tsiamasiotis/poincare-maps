use equilibrium::{Harmonic, HarmonicCache};
use rsl_interpolation::Accelerator;
use safe_unwrap::safe_unwrap;
use utils::{py_eval_harmonic, py_get_numpy1D, py_repr_impl};
use utils::{py_get_float, py_get_path, py_get_typ};

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::PyEqError;

#[derive(Clone)]
#[pyclass(name = "Harmonic")]
pub struct PyHarmonic(pub Harmonic);

#[pymethods]
impl PyHarmonic {
    /// Creates a new PyHarmonic wrapper object.
    #[new]
    pub fn new(path: &str, typ: &str, m: f64, n: f64, phase: f64) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Harmonic::from_dataset(&path, typ, m, n, phase)?))
    }
}

impl From<&Harmonic> for PyHarmonic {
    fn from(harmonic: &Harmonic) -> Self {
        safe_unwrap!(
            "If harmonic exists, Pyharmonic::new() cannot fail",
            PyHarmonic::new(
                safe_unwrap!("file already opened", harmonic.path.to_str()),
                harmonic.typ.as_str(),
                harmonic.m,
                harmonic.n,
                harmonic.phase,
            )
        )
    }
}

py_repr_impl!(PyHarmonic);
py_get_typ!(PyHarmonic);
py_get_path!(PyHarmonic);
py_get_float!(PyHarmonic, psip_wall);
py_get_float!(PyHarmonic, m);
py_get_float!(PyHarmonic, n);
py_get_float!(PyHarmonic, phase);
py_get_numpy1D!(PyHarmonic, psip_data);
py_get_numpy1D!(PyHarmonic, a_data);
py_eval_harmonic!(PyHarmonic, h);
py_eval_harmonic!(PyHarmonic, dh_dpsip);
py_eval_harmonic!(PyHarmonic, dh_dtheta);
py_eval_harmonic!(PyHarmonic, dh_dzeta);
py_eval_harmonic!(PyHarmonic, dh_dt);
