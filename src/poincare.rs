use crate::Particle;
use crate::{Bfield, Current, Qfactor};

use ndarray::{Array1, Array2, Axis};
use numpy::{PyArray2, ToPyArray};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pyclass]
pub struct Poincare {
    pub particles: Vec<Particle>,
    pub angles: Array2<f64>,
    pub fluxes: Array2<f64>,
}

#[pymethods]
impl Poincare {
    #[new]
    pub fn new() -> Self {
        Poincare {
            particles: Vec::new(),
            angles: Array2::zeros((1, 1)),
            fluxes: Array2::zeros((1, 1)),
        }
    }

    #[pyo3(name = "run")]
    fn run_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        angle: &str,
        intersection: f64,
        turns: usize,
    ) -> PyResult<()> {
        dbg!(self.particles.len());
        self.angles = Array2::zeros((self.particles.len(), turns));
        self.fluxes = Array2::zeros((self.particles.len(), turns));

        match self.run(qfactor, bfield, current, angle, intersection, turns) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    pub fn add_particle(&mut self, particle: &Particle) {
        self.particles.push(particle.to_owned())
    }

    pub fn get_angles<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.angles.to_pyarray(py)
    }

    pub fn get_fluxes<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.fluxes.to_pyarray(py)
    }
}

impl Poincare {
    fn run(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        angle: &str,
        intersection: f64,
        turns: usize,
    ) -> PyResult<()> {
        for p in self.particles.iter_mut() {
            match p.run_henon_py(qfactor, bfield, current, angle, intersection, turns) {
                Ok(()) => (),
                Err(err) => return Err(PyTypeError::new_err(err.to_string())),
            };
            dbg!(self.angles.view());
            dbg!(&p.theta.len());
            self.angles
                .push(Axis(1), Array1::from_vec(p.theta.to_owned()).view())
                .unwrap(); // FIXME:
            self.fluxes
                .push(Axis(1), Array1::from_vec(p.psi.to_owned()).view())
                .unwrap(); // FIXME:
        }
        Ok(())
    }
}
