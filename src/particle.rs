use std::time::{Duration, Instant};

use crate::solver::Solver;
use crate::Result;
use crate::{Bfield, Current, Qfactor};
use crate::{InitialConditions, State};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Initial capacity of the vectors that store the evolution of the particle.
const VEC_INIT_CAPACITY: usize = 2000;

/// The initial time step for the RFK45 adaptive step method. Should be small
/// enough to account for fast particles. The value is empirical.
const RKF45_FIRST_STEP: f64 = 1e-4;

#[pyclass]
pub struct Particle {
    /// The initial (θ, ψ_p, ρ, ζ, μ) of the particle.
    initial: InitialConditions,
    /// The current state of the particle.
    state: State,
    /// The current RK45 state.
    solver: Solver,
    /// The calculated evaluation times.
    #[pyo3(get)]
    t: Vec<f64>,
    /// The calculated θ values.
    #[pyo3(get)]
    theta: Vec<f64>,
    /// The calculated ψ_p values.
    #[pyo3(get)]
    psip: Vec<f64>,
    /// The calculated ρ values.
    #[pyo3(get)]
    rho: Vec<f64>,
    /// The calculated ζ values.
    #[pyo3(get)]
    zeta: Vec<f64>,
    /// The calculated Pζ values.
    #[pyo3(get)]
    pzeta: Vec<f64>,
    /// The calculated Pθ values.
    #[pyo3(get)]
    ptheta: Vec<f64>,
    /// The calculated ψ values.
    #[pyo3(get)]
    psi: Vec<f64>,
    initial_energy: f64,
    final_energy: f64,
    calculation_time: Duration,
}

#[pymethods]
impl Particle {
    #[new]
    pub fn new(initial: InitialConditions) -> Self {
        Self {
            initial: initial.to_owned(),
            state: State::new_init(&initial),
            solver: Solver::default(),
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

    #[pyo3(name = "run_henon_zeta")]
    pub fn run_henon_zeta_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        angle: f64,
        turns: usize,
    ) -> PyResult<()> {
        match self.run_henon_zeta(qfactor, bfield, current, angle, turns) {
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
            self.solver = Solver::default();
            self.solver.init(&self.state);
            self.solver.start(h, qfactor, bfield, current)?;
            h = self.solver.calculate_optimal_step(h);
            self.state = self.solver.next_state(h);

            self.state.evaluate(qfactor, current, bfield)?;
            self.update_vecs();
        }
        self.calculation_time = start.elapsed();
        self.shrink_vecs();

        self.final_energy = self.state.energy();

        Ok(())
    }

    pub fn run_henon_zeta(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        angle: f64,
        turns: usize,
    ) -> Result<()> {
        use std::f64::consts::TAU;
        self.state.evaluate(qfactor, current, bfield)?;
        self.initial_energy = self.state.energy();

        let mut h = RKF45_FIRST_STEP;
        let angle = angle % TAU;

        self.calculation_time = Duration::ZERO;
        let start = Instant::now();

        while self.zeta.len() <= turns {
            self.solver = Solver::default();
            self.solver.init(&self.state);
            self.solver.start(h, qfactor, bfield, current)?;
            h = self.solver.calculate_optimal_step(h);
            self.state = self.solver.next_state(h);

            self.state.evaluate(qfactor, current, bfield)?;

            let zeta_old = self.solver.state1.zeta;
            let zeta_new = self.state.zeta;
            if intersected(zeta_old, zeta_new, angle) {
                // TODO: calculate state with Henon's trick
                self.update_vecs();
            }
        }

        // TODO: tighten the restriction
        assert!(self
            .zeta
            .windows(2)
            .map(|w| (w[1] - w[0]).abs() - TAU <= 1e-2)
            .all(|b| b));

        self.calculation_time = start.elapsed();
        self.shrink_vecs();

        self.final_energy = self.state.energy();

        Ok(())
    }

    fn shrink_vecs(&mut self) {
        self.t.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.psi.shrink_to_fit();
    }

    fn update_vecs(&mut self) {
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

/// Checks when an angle has intersected with the surface at `angle`.
/// (source: seems to work)
fn intersected(old_angle: f64, new_angle: f64, surface_angle: f64) -> bool {
    ((new_angle - surface_angle) / 2.0).sin() * ((old_angle - surface_angle) / 2.0).sin() <= 0.0
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
    RKF45_FIRST_STEP
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let evolution = if self.t.is_empty() {
            "not calculated".to_string()
        } else {
            format!("{} steps", self.t.len())
        };

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
            .field("ψ-acc", &self.state.xacc)
            .field("θ-acc", &self.state.xacc)
            .field("Evolution", &evolution)
            .field("Initial energy", &self.initial_energy)
            .field("Final energy  ", &self.final_energy)
            .field("Time", &self.calculation_time)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::PathBuf;

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
        let mut particle = Particle::new(initial);
        particle
            .run_ode(&qfactor, &bfield, &current, (0.0, 210.0), 50000)
            .unwrap();
        dbg!(&particle);
    }
}
