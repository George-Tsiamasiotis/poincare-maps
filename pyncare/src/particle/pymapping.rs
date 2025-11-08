use particle::MappingParameters;
use utils::{py_debug_impl, py_repr_impl};

use pyo3::prelude::*;

#[pyclass(name = "MappingParameters")]
pub struct PyMappingParameters(pub MappingParameters);

#[pymethods]
impl PyMappingParameters {
    #[new]
    pub fn new(section: &str, alpha: f64, intersections: usize) -> Self {
        let section = match section {
            "theta" => particle::PoincareSection::ConstTheta,
            "zeta" => particle::PoincareSection::ConstZeta,
            _ => panic!("mapping angle must be 'theta' or 'zeta'"),
        };
        Self(MappingParameters::new(section, alpha, intersections))
    }

    #[getter]
    pub fn get_section(&self) -> String {
        format!("{:?}", self.0.section)
    }

    #[getter]
    pub fn get_alpha(&self) -> f64 {
        self.0.alpha
    }

    #[getter]
    pub fn get_intersections(&self) -> usize {
        self.0.intersections
    }
}

py_debug_impl!(PyMappingParameters);
py_repr_impl!(PyMappingParameters);
