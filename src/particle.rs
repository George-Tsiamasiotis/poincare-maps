use std::time::{Duration, Instant};

use crate::solver::Solver;
use crate::Result;
use crate::{Bfield, Current, Qfactor};
use crate::{InitialConditions, State};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Initial capacity of the vectors that store the evolution of the particle.
const VEC_INIT_CAPACITY: usize = 2000;

/// The initial time step for the RKF45 adaptive step method. Should be small
/// enough to account for fast particles. The value is empirical.
const RKF45_FIRST_STEP: f64 = 1e-4;

#[pyclass]
pub struct Particle {
    /// The initial (θ, ψ_p, ρ, ζ, μ) of the particle.
    initial: InitialConditions,
    /// The current state of the particle.
    state: State,
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
            let mut solver = Solver::default();
            solver.init(&self.state);
            solver.start(h, qfactor, bfield, current)?;
            h = solver.calculate_optimal_step(h);
            self.state = solver.next_state(h);

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

        self.calculation_time = Duration::ZERO;
        let start = Instant::now();

        while self.zeta.len() <= turns {
            let mut solver = Solver::default();
            solver.init(&self.state);
            solver.start(h, qfactor, bfield, current)?;
            h = solver.calculate_optimal_step(h);

            // State before the intersection
            let old_state = self.state.clone();

            let mut next_state = solver.next_state(h);
            next_state.evaluate(qfactor, current, bfield)?;

            if intersected(old_state.zeta, next_state.zeta, angle) {
                // Setup new system (6) and (9)
                let kappa = 1.0 / old_state.zeta_dot;
                let dtheta_dzeta = kappa * old_state.theta_dot;
                let dpsip_dzeta = kappa * old_state.psip_dot;
                let drho_dzeta = kappa * old_state.rho_dot;
                let dt_dzeta = kappa;

                let mod_state = State {
                    t: old_state.zeta,
                    zeta: old_state.t,
                    theta_dot: dtheta_dzeta,
                    psip_dot: dpsip_dzeta,
                    rho_dot: drho_dzeta,
                    zeta_dot: dt_dzeta,
                    ..old_state
                };

                let direction = (next_state.t - old_state.t).signum();
                // NOTE: unsure about where tho MOD should happen
                let dz = direction * (angle - old_state.zeta % TAU);

                let mut solver = Solver::default();
                solver.init(&mod_state);
                solver.start(dz, qfactor, bfield, current)?;
                let new_mod_state = solver.next_state(dz);

                let kappa = 1.0;
                let dtheta_dt = kappa * new_mod_state.theta_dot;
                let dpsip_dt = kappa * new_mod_state.psip_dot;
                let drho_dt = kappa * new_mod_state.rho_dot;
                let dt_dt = kappa;
                let intersection_state = State {
                    t: new_mod_state.zeta,
                    zeta: new_mod_state.t,
                    theta_dot: dtheta_dt,
                    psip_dot: dpsip_dt,
                    rho_dot: drho_dt,
                    zeta_dot: dt_dt,
                    ..new_mod_state
                };

                self.state = intersection_state.clone();
                self.state.evaluate(qfactor, current, bfield)?;
                self.update_vecs();

                let temp_step = next_state.t - intersection_state.t;
                let mut solver = Solver::default();
                solver.init(&self.state);
                solver.start(temp_step, qfactor, bfield, current)?;
                self.state = solver.next_state(temp_step);
                self.state.evaluate(qfactor, current, bfield)?;
                h = RKF45_FIRST_STEP;
            } else {
                self.state = next_state;
            }
        }

        // TODO: tighten the restriction
        assert!(self
            .zeta
            .windows(2)
            .map(|w| (w[1] - w[0]).abs() - TAU <= 1e-9)
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
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    let cond = ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() < 0.0;
    if cond {
        assert!(diff1.sin() >= 0.0);
        assert!(diff2.sin() <= 0.0);
    }
    cond
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
        let mut particle = Particle::new(initial);
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
        let mut particle = Particle::new(initial);
        particle
            .run_henon_zeta(&qfactor, &bfield, &current, PI / 2.0, 10)
            .unwrap();
        dbg!(&particle);
    }
}
