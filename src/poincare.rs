use std::time::Duration;

use crate::MapError;
use crate::Particle;
use crate::PoincareParameters;
use crate::Result;
use crate::Surface;
use crate::solver::henon;
use crate::{__repr__, numpy_getter_2D};
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
    pub params: PoincareParameters,
    #[pyo3(get)]
    pub timed_out_particles: usize,
    #[pyo3(get)]
    pub completed_particles: usize,
    #[pyo3(get)]
    pub escpaped_particles: usize,
    #[pyo3(get)]
    pub max_steps: usize,
    #[pyo3(get)]
    pub max_duration: Duration,
}

#[pymethods]
impl Poincare {
    #[new]
    pub fn new(params: PoincareParameters) -> Self {
        Poincare {
            particles: Vec::new(),
            angles: Array2::zeros((1, 1)),
            fluxes: Array2::zeros((1, 1)),
            params,
            timed_out_particles: 0,
            completed_particles: 0,
            escpaped_particles: 0,
            max_steps: 0,
            max_duration: Duration::default(),
        }
    }

    #[pyo3(name = "run")]
    #[coverage(off)]
    pub fn run_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> PyResult<()> {
        match self.run(qfactor, bfield, current, per) {
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
}

__repr__!(Poincare);
numpy_getter_2D!(Poincare, angles);
numpy_getter_2D!(Poincare, fluxes);

impl Poincare {
    #[allow(clippy::too_many_arguments)]
    pub fn run(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> Result<()> {
        self.angles = Array2::zeros((0, self.params.turns));
        self.fluxes = Array2::zeros((0, self.params.turns));

        let style = ProgressStyle::with_template(
            "[{elapsed_precise}] {wide_bar:.cyan/blue} {spinner} {pos:>4}/{len:4} {msg}",
        )
        .unwrap();
        let pbar = ProgressBar::new(self.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        pbar.force_draw();

        // Start a new thread for each particle
        self.particles.par_iter_mut().try_for_each(|p| {
            {
                match henon::run_henon(
                    &mut p.evolution,
                    qfactor,
                    bfield,
                    current,
                    per,
                    &self.params,
                ) {
                    // Discard particles that hit the wall instead of panicking.
                    Err(err) if matches!(err, MapError::DomainError(..)) => {
                        // p.status = IntegrationStatus::Escaped;
                        Ok(())
                    }
                    // Discard particles that took too long to integrate
                    Err(err) if matches!(err, MapError::OrbitTimeout(..)) => {
                        // p.status = IntegrationStatus::TimedOut;
                        Ok(())
                    }
                    Err(err) => Err(err),
                    Ok(_) => Ok(()),
                }
            }
            .inspect(|()| pbar.inc(1))
        })?;

        // Store points
        for p in self.particles.iter_mut() {
            // Do not store particles that didn't complete the integration.
            if p.evolution.time.len() != self.params.turns {
                continue;
            }
            match self.params.surface {
                Surface::ConstZeta => {
                    self.angles
                        .push_row(Array1::from_vec(p.evolution.zeta.clone()).view())
                        .unwrap();
                    self.fluxes
                        .push_row(Array1::from_vec(p.evolution.psip.clone()).view())
                        .unwrap()
                }
                Surface::ConstTheta => {
                    self.angles
                        .push_row(Array1::from_vec(p.evolution.theta.clone()).view())
                        .unwrap();
                    self.fluxes
                        .push_row(Array1::from_vec(p.evolution.psi.clone()).view())
                        .unwrap()
                }
            }
        }

        pbar.finish_with_message("Done");
        self.statistics();
        Ok(())
    }

    fn statistics(&mut self) {
        // for p in self.particles.iter() {
        //     match p.status {
        //         IntegrationStatus::Initialized => (),
        //         IntegrationStatus::Integrated => self.completed_particles += 1,
        //         IntegrationStatus::Escaped => self.escpaped_particles += 1,
        //         IntegrationStatus::TimedOut => self.timed_out_particles += 1,
        //     }
        //     self.max_steps = p.steps_taken.max(self.max_steps);
        //     self.max_duration = p.calculation_time.max(self.max_duration);
        // }
    }

    /// Adds a particle to be calculated.
    pub fn add_particle(&mut self, particle: &Particle) {
        self.particles.push(particle.to_owned())
    }

    /// Returns the stored particles.
    pub fn get_particles(&self) -> Vec<Particle> {
        self.particles.clone()
    }
}

impl std::fmt::Debug for Poincare {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Poincare")
            .field("number of particles", &self.particles.len())
            .field("params", &self.params)
            .field("timed out particles", &self.timed_out_particles)
            .field("completed particles", &self.completed_particles)
            .field("escpaped_particles", &self.escpaped_particles)
            .field("max duration", &self.max_duration)
            .field("max steps", &self.max_steps)
            .finish()
    }
}
