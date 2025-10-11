use crate::Particle;
use crate::{Bfield, Current, Qfactor};

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array1, Array2};
use numpy::{PyArray2, ToPyArray};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pyclass]
pub struct Poincare {
    pub particles: Vec<Particle>,
    pub angles: Array2<f64>,
    pub fluxes: Array2<f64>,
    #[pyo3(get)]
    pub angle: String,
    #[pyo3(get)]
    pub intersection: f64,
}

#[pymethods]
impl Poincare {
    #[new]
    pub fn new() -> Self {
        Poincare {
            particles: Vec::new(),
            angles: Array2::zeros((1, 1)),
            fluxes: Array2::zeros((1, 1)),
            angle: "".into(),
            intersection: f64::NAN,
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
        self.angles = Array2::zeros((0, turns));
        self.fluxes = Array2::zeros((0, turns));
        self.angle = angle.into();
        self.intersection = intersection;

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
        let style = ProgressStyle::with_template(
            "[{elapsed_precise}] {wide_bar:.cyan/blue} {spinner} {pos:>4}/{len:4} {msg}",
        )
        .unwrap();
        let pbar = ProgressBar::new(self.particles.len() as u64).with_style(style);
        for p in self.particles.iter_mut() {
            match p.run_henon_py(qfactor, bfield, current, angle, intersection, turns) {
                Ok(()) => (),
                Err(err) => return Err(PyTypeError::new_err(err.to_string())),
            };
            match angle {
                "theta" => {
                    self.angles
                        .push_row(Array1::from_vec(p.zeta.to_owned()).view())
                        .unwrap();
                    self.fluxes
                        .push_row(Array1::from_vec(p.psip.to_owned()).view())
                        .unwrap();
                }
                "zeta" => {
                    self.angles
                        .push_row(Array1::from_vec(p.theta.to_owned()).view())
                        .unwrap();
                    self.fluxes
                        .push_row(Array1::from_vec(p.psi.to_owned()).view())
                        .unwrap();
                }
                _ => unreachable!(),
            }
            pbar.inc(1);
        }
        pbar.finish_with_message("Done");

        Ok(())
    }
}
