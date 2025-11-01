use crate::repr_impl;
use particle::Mapping;

use pyo3::prelude::*;

#[pyclass(name = "Mapping")]
pub struct PyMapping {
    pub mapping: Mapping,
}

#[pymethods]
impl PyMapping {
    #[new]
    pub fn new(section: &str, alpha: f64, intersections: usize) -> Self {
        let section = match section {
            "theta" => particle::PoincareSection::ConstTheta,
            "zeta" => particle::PoincareSection::ConstZeta,
            _ => panic!("mapping angle must be 'theta' or 'zeta'"),
        };
        Self {
            mapping: Mapping::new(section, alpha, intersections),
        }
    }

    #[getter]
    pub fn get_section(&self) -> String {
        format!("{:?}", self.mapping.section)
    }

    #[getter]
    pub fn get_alpha(&self) -> f64 {
        self.mapping.alpha
    }

    #[getter]
    pub fn get_intersections(&self) -> usize {
        self.mapping.intersections
    }
}

repr_impl!(PyMapping);

impl std::fmt::Debug for PyMapping {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.mapping.fmt(f)
    }
}
