use std::time::Duration;

use crate::Particle;
use crate::Result;
use crate::solver::henon;
use crate::{Bfield, Current, Perturbation, Qfactor};

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array1, Array2};
use numpy::{PyArray2, ToPyArray};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

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
        #[cfg(feature = "rk45")]
        compile_error!("Feature rk45 must be disabled for Poincare maps.");

        Poincare {
            particles: Vec::new(),
            angles: Array2::zeros((1, 1)),
            fluxes: Array2::zeros((1, 1)),
            angle: "".into(),
            intersection: f64::NAN,
        }
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(name = "run")]
    #[coverage(off)]
    pub fn run_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        angle: &str,
        intersection: f64,
        turns: usize,
    ) -> PyResult<()> {
        match self.run(qfactor, bfield, current, per, angle, intersection, turns) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    /// Adds a particle to be calculated.
    #[pyo3(name = "get_particles")]
    #[coverage(off)]
    pub fn get_particles_py(&self) -> Vec<Particle> {
        self.particles.clone()
    }

    /// Returns the stored particles in a list.
    #[pyo3(name = "add_particle")]
    #[coverage(off)]
    pub fn add_particle_py(&mut self, particle: &Particle) {
        self.add_particle(particle);
    }

    /// Returns the calculated angles as a 2D numpy array, 1 row corresponding to 1 particle.
    #[pyo3(name = "get_angles")]
    #[coverage(off)]
    pub fn get_angles_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.get_angles().to_pyarray(py)
    }

    /// Returns the calculated fluxes as a 2D numpy array, 1 row corresponding to 1 particle.
    #[pyo3(name = "get_fluxes")]
    #[coverage(off)]
    pub fn get_fluxes_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.fluxes.to_pyarray(py)
    }
}

impl Poincare {
    #[allow(clippy::too_many_arguments)]
    pub fn run(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        angle: &str,
        intersection: f64,
        turns: usize,
    ) -> Result<()> {
        self.angles = Array2::zeros((0, turns));
        self.fluxes = Array2::zeros((0, turns));
        self.angle = angle.into();
        self.intersection = intersection;

        let style = ProgressStyle::with_template(
            "[{elapsed_precise}] {wide_bar:.cyan/blue} {spinner} {pos:>4}/{len:4} {msg}",
        )
        .unwrap();
        let pbar = ProgressBar::new(self.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        pbar.force_draw();

        // Start a new thread for each particle
        self.particles.par_iter_mut().try_for_each(|p| {
            henon::run_henon(p, qfactor, bfield, current, per, angle, intersection, turns)
                .inspect(|()| pbar.inc(1))
        })?;

        // Store points
        for p in self.particles.iter_mut() {
            match angle {
                "theta" => {
                    self.angles
                        .push_row(Array1::from_vec(p.zeta.clone()).view())
                        .unwrap();
                    self.fluxes
                        .push_row(Array1::from_vec(p.psip.clone()).view())
                        .unwrap()
                }
                "zeta" => {
                    self.angles
                        .push_row(Array1::from_vec(p.theta.clone()).view())
                        .unwrap();
                    self.fluxes
                        .push_row(Array1::from_vec(p.psi.clone()).view())
                        .unwrap()
                }
                _ => unreachable!(),
            }
        }

        pbar.finish_with_message("Done");

        Ok(())
    }

    /// Adds a particle to be calculated.
    pub fn add_particle(&mut self, particle: &Particle) {
        self.particles.push(particle.to_owned())
    }

    /// Returns the stored particles.
    pub fn get_particles(&self) -> Vec<Particle> {
        self.particles.clone()
    }

    /// Returns the calculated angles as a 2D array, 1 row corresponding to 1 particle.
    pub fn get_angles(&self) -> Array2<f64> {
        self.angles.clone()
    }

    /// Returns the calculated fluxes as a 2D array, 1 row corresponding to 1 particle.
    pub fn get_fluxes(&self) -> Array2<f64> {
        self.angles.clone()
    }
}

impl Default for Poincare {
    fn default() -> Self {
        Self::new()
    }
}
