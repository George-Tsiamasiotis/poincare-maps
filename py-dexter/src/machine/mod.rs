//! Defines `dexter-machine` helper objects and exports machine objects and helper types

mod bfield;
mod current;
mod geometry;
mod mode;
mod perturbation;
mod qfactor;

pub use bfield::*;
pub use current::*;
pub use geometry::*;
pub use mode::*;
pub use perturbation::*;
pub use qfactor::*;

use crate::*;
use dexter::dexter_machine::*;
use pyo3::{prelude::*, types::PyType};

#[pyclass(name = "_PyLastClosedFluxSurface", frozen, immutable_type)]
pub struct PyLastClosedFluxSurface(pub LastClosedFluxSurface);

#[pymethods]
impl PyLastClosedFluxSurface {
    #[classmethod]
    pub fn toroidal(_: &Bound<'_, PyType>, value: f64) -> PyResult<Self> {
        Ok(Self(LastClosedFluxSurface::Toroidal(value)))
    }

    #[classmethod]
    pub fn poloidal(_: &Bound<'_, PyType>, value: f64) -> PyResult<Self> {
        Ok(Self(LastClosedFluxSurface::Poloidal(value)))
    }

    #[getter]
    pub fn value(&self) -> f64 {
        self.0.value()
    }

    #[getter]
    pub fn kind(&self) -> String {
        match self.0 {
            LastClosedFluxSurface::Toroidal(_) => "Toroidal".into(),
            LastClosedFluxSurface::Poloidal(_) => "Poloidal".into(),
        }
    }
}

impl_py_repr!(PyLastClosedFluxSurface, simple);
