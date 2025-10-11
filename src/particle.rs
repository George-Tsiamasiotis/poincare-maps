use std::time::{Duration, Instant};

use crate::Result;
use crate::solver::Solver;
use crate::solver::henon;
use crate::{Bfield, Current, Qfactor};
use crate::{InitialConditions, State};

use numpy::{PyArray1, ToPyArray};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Initial capacity of the vectors that store the evolution of the particle.
const VEC_INIT_CAPACITY: usize = 2000;

#[pyclass]
#[derive(Clone)]
pub struct Particle {
    /// The initial (θ, ψ_p, ρ, ζ, μ) of the particle.
    #[pyo3(get)]
    pub initial: InitialConditions,
    /// The current state of the particle.
    pub(crate) state: State,
    /// The calculated evaluation times.
    #[pyo3(get)]
    pub t: Vec<f64>,
    /// The calculated θ values.
    #[pyo3(get)]
    pub theta: Vec<f64>,
    /// The calculated ψ_p values.
    #[pyo3(get)]
    pub psip: Vec<f64>,
    /// The calculated ρ values.
    #[pyo3(get)]
    pub rho: Vec<f64>,
    /// The calculated ζ values.
    #[pyo3(get)]
    pub zeta: Vec<f64>,
    /// The calculated Pζ values.
    #[pyo3(get)]
    pub pzeta: Vec<f64>,
    /// The calculated Pθ values.
    #[pyo3(get)]
    pub ptheta: Vec<f64>,
    /// The calculated ψ values.
    #[pyo3(get)]
    pub psi: Vec<f64>,
    /// The starting energy of the particle.
    #[pyo3(get)]
    pub initial_energy: f64,
    /// The final energy of the particle.
    #[pyo3(get)]
    pub final_energy: f64,
    /// The total compute time.
    pub(crate) calculation_time: Duration,
    /// The total number of integration steps.
    pub(crate) steps_taken: usize,
}

#[pymethods]
impl Particle {
    /// Creates a new particle from the initial conditions.
    #[new]
    pub fn new(initial: &InitialConditions) -> Self {
        let mut t = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut theta = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut psip = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut rho = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut zeta = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut psi = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut ptheta = Vec::with_capacity(VEC_INIT_CAPACITY);
        let mut pzeta = Vec::with_capacity(VEC_INIT_CAPACITY);
        t.push(initial.t0);
        theta.push(initial.theta0);
        psip.push(initial.psip0);
        rho.push(initial.rho0);
        zeta.push(initial.zeta0);
        psi.push(f64::NAN);
        ptheta.push(f64::NAN);
        pzeta.push(f64::NAN);
        Self {
            initial: initial.to_owned(),
            state: State::new_init(&initial),
            t,
            theta,
            psip,
            rho,
            zeta,
            psi,
            ptheta,
            pzeta,
            initial_energy: f64::NAN,
            final_energy: f64::NAN,
            calculation_time: Duration::ZERO,
            steps_taken: 0,
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }

    /// Calculates the particle's trajectory.
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
    /// Calculates the particle's trajectory, by also stepping exactly at the `intersection`
    /// surface, with respect to `angle`, which can be either "theta" or "zeta".
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

    #[getter]
    pub fn t<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.t.to_pyarray(py)
    }

    #[getter]
    pub fn theta<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.theta.to_pyarray(py)
    }

    #[getter]
    pub fn psip<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.psip.to_pyarray(py)
    }

    #[getter]
    pub fn rho<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.rho.to_pyarray(py)
    }

    #[getter]
    pub fn zeta<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.zeta.to_pyarray(py)
    }

    #[getter]
    pub fn ptheta<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.ptheta.to_pyarray(py)
    }

    #[getter]
    pub fn pzeta<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.pzeta.to_pyarray(py)
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
