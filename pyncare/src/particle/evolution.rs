use particle::Evolution;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use crate::to_numpy1D_impl;

#[pyclass]
#[derive(Debug, Clone)]
#[pyo3(name = "Evolution")]
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

to_numpy1D_impl!(PyEvolution, evolution, time);
to_numpy1D_impl!(PyEvolution, evolution, theta);
to_numpy1D_impl!(PyEvolution, evolution, psip);
to_numpy1D_impl!(PyEvolution, evolution, rho);
to_numpy1D_impl!(PyEvolution, evolution, zeta);
to_numpy1D_impl!(PyEvolution, evolution, psi);
to_numpy1D_impl!(PyEvolution, evolution, ptheta);
to_numpy1D_impl!(PyEvolution, evolution, pzeta);
