use std::time::{Duration, Instant};

use crate::solver::henon;
use crate::solver::Solver;
use crate::Result;
use crate::{Bfield, Current, Qfactor};
use crate::{InitialConditions, State};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Initial capacity of the vectors that store the evolution of the particle.
const VEC_INIT_CAPACITY: usize = 2000;

#[pyclass]
pub struct Particle {
    /// The initial (θ, ψ_p, ρ, ζ, μ) of the particle.
    pub(crate) initial: InitialConditions,
    /// The current state of the particle.
    pub(crate) state: State,
    /// The calculated evaluation times.
    #[pyo3(get)]
    pub(crate) t: Vec<f64>,
    /// The calculated θ values.
    #[pyo3(get)]
    pub(crate) theta: Vec<f64>,
    /// The calculated ψ_p values.
    #[pyo3(get)]
    pub(crate) psip: Vec<f64>,
    /// The calculated ρ values.
    #[pyo3(get)]
    pub(crate) rho: Vec<f64>,
    /// The calculated ζ values.
    #[pyo3(get)]
    pub(crate) zeta: Vec<f64>,
    /// The calculated Pζ values.
    #[pyo3(get)]
    pub(crate) pzeta: Vec<f64>,
    /// The calculated Pθ values.
    #[pyo3(get)]
    pub(crate) ptheta: Vec<f64>,
    /// The calculated ψ values.
    #[pyo3(get)]
    pub(crate) psi: Vec<f64>,
    pub(crate) initial_energy: f64,
    pub(crate) final_energy: f64,
    pub(crate) calculation_time: Duration,
    /// The total number of steps taken.
    pub(crate) steps_taken: usize,
}

#[pymethods]
impl Particle {
    #[new]
    pub fn new(initial: &InitialConditions) -> Self {
        Self {
            initial: initial.to_owned(),
            state: State::new_init(&initial),
            t: Vec::with_capacity(VEC_INIT_CAPACITY),
            theta: Vec::with_capacity(VEC_INIT_CAPACITY),
            psip: Vec::with_capacity(VEC_INIT_CAPACITY),
            rho: Vec::with_capacity(VEC_INIT_CAPACITY),
            zeta: Vec::with_capacity(VEC_INIT_CAPACITY),
            pzeta: Vec::with_capacity(VEC_INIT_CAPACITY),
            ptheta: Vec::with_capacity(VEC_INIT_CAPACITY),
            psi: Vec::with_capacity(VEC_INIT_CAPACITY),
            initial_energy: f64::NAN,
            final_energy: f64::NAN,
            calculation_time: Duration::ZERO,
            steps_taken: 0,
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }

    #[pyo3(name = "run_ode")]
    pub fn run_ode_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        t_eval: (f64, f64),
        steps: usize,
    ) -> PyResult<()> {
        match self.run_ode(qfactor, bfield, current, t_eval, steps) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    #[pyo3(name = "run_henon")]
    pub fn run_henon_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        angle: &str,
        intersection: f64,
        turns: usize,
    ) -> PyResult<()> {
        match henon::run_henon(self, qfactor, bfield, current, angle, intersection, turns) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }
}

impl Particle {
    pub fn run_ode(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        t_eval: (f64, f64),
        steps: usize,
    ) -> Result<()> {
        self.state.evaluate(qfactor, current, bfield)?;
        self.initial_energy = self.state.energy();

        let mut h = first_step(t_eval, steps);

        self.calculation_time = Duration::ZERO;
        let start = Instant::now();

        while self.state.t < t_eval.1 {
            let mut solver = Solver::default();
            solver.init(&self.state);
            solver.start(h, qfactor, bfield, current)?;
            h = solver.calculate_optimal_step(h);
            self.state = solver.next_state(h);

            self.state.evaluate(qfactor, current, bfield)?;
            self.update_vecs();
            self.steps_taken += 1;
        }
        self.calculation_time = start.elapsed();
        self.shrink_vecs();

        self.final_energy = self.state.energy();

        Ok(())
    }

    pub(crate) fn shrink_vecs(&mut self) {
        self.t.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.psi.shrink_to_fit();
    }

    pub(crate) fn update_vecs(&mut self) {
        self.t.push(self.state.t);
        self.theta.push(self.state.theta);
        self.psip.push(self.state.psip);
        self.rho.push(self.state.rho);
        self.zeta.push(self.state.zeta);
        self.pzeta.push(self.state.pzeta);
        self.ptheta.push(self.state.ptheta);
        self.psi.push(self.state.psi);
    }
}

#[cfg(feature = "rk45")]
/// Returns the standard time step, which will be constant throughout the integration.
fn first_step(t_eval: (f64, f64), steps: usize) -> f64 {
    let t0 = t_eval.0;
    let tf = t_eval.1;
    let step = (tf - t0) / (steps as f64);

    step
}

#[cfg(not(feature = "rk45"))]
#[allow(unused_variables)]
fn first_step(t_eval: (f64, f64), steps: usize) -> f64 {
    crate::solver::RKF45_FIRST_STEP
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("Initial", &self.initial)
            .field(
                "Interval",
                &format!(
                    "[{:.5}, {:.5}]",
                    self.initial.t0,
                    self.t.last().copied().unwrap_or(f64::NAN)
                ),
            )
            .field("Parallel energy", &self.state.parallel_energy())
            .field("Perpendicular energy", &self.state.perpendicular_energy())
            .field("ψ-acc", &self.state.xacc)
            .field("θ-acc", &self.state.xacc)
            .field("Steps taken", &self.steps_taken)
            .field("Steps stored", &self.t.len())
            .field("Initial energy", &self.initial_energy)
            .field("Final energy  ", &self.final_energy)
            .field("Time", &self.calculation_time)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::f64::consts::TAU;
    use std::{f64::consts::PI, path::PathBuf};

    #[test]
    fn test_particle_run_ode() {
        let path = PathBuf::from("./data.nc");
        let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
        let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
        let current = Current::from_dataset(&path, "akima").unwrap();

        // Normal passing particle
        let initial = InitialConditions {
            t0: 0.0,
            theta0: 3.14,
            psip0: 0.05,
            rho0: 0.05,
            zeta0: 0.1,
            mu: 1e-4,
        };
        let mut particle = Particle::new(&initial);
        particle
            .run_ode(&qfactor, &bfield, &current, (0.0, 210.0), 50000)
            .unwrap();
        dbg!(&particle);
    }

    #[test]
    fn test_particle_run_henon_zeta() {
        let path = PathBuf::from("./data.nc");
        let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
        let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
        let current = Current::from_dataset(&path, "akima").unwrap();

        // Normal passing particle
        let initial = InitialConditions {
            t0: 0.0,
            theta0: 3.14,
            psip0: 0.02,
            rho0: 0.05,
            zeta0: 0.0,
            mu: 1e-4,
        };
        let mut particle = Particle::new(&initial);
        particle
            .run_henon_py(&qfactor, &bfield, &current, "zeta", PI / 2.0, 10)
            .unwrap();
        dbg!(&particle.zeta.iter().map(|z| z % TAU).collect::<Vec<f64>>());
        dbg!(&particle);

        let mut particle = Particle::new(&initial);
        particle
            .run_henon_py(&qfactor, &bfield, &current, "theta", PI / 2.0, 10)
            .unwrap();
        dbg!(&particle.zeta.iter().map(|z| z % TAU).collect::<Vec<f64>>());
        dbg!(&particle);
    }
}
