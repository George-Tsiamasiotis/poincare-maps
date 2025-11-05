use particle::Evolution;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use utils::{repr_impl, to_numpy1D_impl};

#[derive(Clone)]
#[pyclass(name = "Evolution")]
pub struct PyEvolution {
    pub evolution: Evolution,
}

impl PyEvolution {
    pub fn from_evolution(evolution: &Evolution) -> Self {
        Self {
            evolution: evolution.clone(),
        }
    }
}

#[pymethods]
impl PyEvolution {
    #[getter]
    pub fn get_duration(&self) -> String {
        format!("{:?}", self.evolution.duration)
    }
}

repr_impl!(PyEvolution);
to_numpy1D_impl!(PyEvolution, evolution, time);
to_numpy1D_impl!(PyEvolution, evolution, theta);
to_numpy1D_impl!(PyEvolution, evolution, psip);
to_numpy1D_impl!(PyEvolution, evolution, rho);
to_numpy1D_impl!(PyEvolution, evolution, zeta);
to_numpy1D_impl!(PyEvolution, evolution, psi);
to_numpy1D_impl!(PyEvolution, evolution, ptheta);
to_numpy1D_impl!(PyEvolution, evolution, pzeta);

impl std::fmt::Debug for PyEvolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.evolution.fmt(f)
    }
}
