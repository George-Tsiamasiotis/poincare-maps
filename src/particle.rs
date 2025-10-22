use crate::__repr__;
use crate::Evolution;
use crate::IntegrationStatus;
use crate::PoincareParameters;
use crate::Result;
use crate::State;
use crate::solver::Solver;
use crate::solver::henon;
use crate::{Bfield, Current, Perturbation, Qfactor};

use crate::consts::RKF45_FIRST_STEP;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Representation of a particle or a magnetic field line, for μ=0 and ρ->0.
#[pyclass]
#[derive(Default, Clone)]
pub struct Particle {
    /// The initial state of the particle.
    pub initial_state: State,
    /// The final state of the particle after the integration.
    pub final_state: State,
    /// The evolution time series of the particle.
    pub evolution: Evolution,
    /// The starting energy of the particle.
    pub initial_energy: f64,
    /// The final energy of the particle.
    pub final_energy: f64,
    /// Status about the particle's integration.
    pub status: IntegrationStatus,
}

#[pymethods]
impl Particle {
    /// Creates a new particle from the initial conditions.
    #[new]
    pub fn new(t0: f64, theta0: f64, psip0: f64, rho0: f64, zeta0: f64, mu: f64) -> Self {
        Self {
            initial_state: State::from_initial(t0, theta0, psip0, rho0, zeta0, mu),
            ..Default::default()
        }
    }

    /// Calculates the particle's trajectory.
    #[pyo3(name = "run_ode")]
    #[coverage(off)]
    pub fn run_ode_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        t_eval: (f64, f64),
    ) -> PyResult<()> {
        match self.run_ode(qfactor, bfield, current, per, t_eval) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    /// Calculates the particle's trajectory, by also stepping exactly at the `intersection`
    /// surface, with respect to `angle`, which can be either "theta" or "zeta".
    #[pyo3(name = "run_henon")]
    #[coverage(off)]
    pub fn run_henon_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        params: &PoincareParameters,
    ) -> PyResult<()> {
        match henon::run_henon(&mut self.evolution, qfactor, bfield, current, per, params) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }
}

__repr__!(Particle);

impl Particle {
    pub fn run_ode(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        t_eval: (f64, f64),
    ) -> Result<()> {
        let initial_state = self
            .initial_state
            .clone()
            .into_evaluated(qfactor, current, bfield, per)?;

        // Tracks the state of the particle at each step.
        let mut state = initial_state.clone();
        let mut h = RKF45_FIRST_STEP;

        while state.time < t_eval.1 {
            // Calculate next state.
            let mut solver = Solver::default();
            solver.init(&state);
            solver.start(h, qfactor, bfield, current, per)?;
            h = solver.calculate_optimal_step(h);
            state = solver.next_state(h);
            state.evaluate(qfactor, current, bfield, per)?;

            self.evolution.push_point(state.as_point());
        }

        // Finalize
        self.final_state = state;
        self.evolution.shrink_to_fit();
        self.status = IntegrationStatus::Integrated;
        self.final_energy = self.final_state.energy();
        Ok(())
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("Initial", &self.initial_state)
            .field(
                "Interval",
                &format!(
                    "[{:.5}, {:.5}]",
                    self.initial_state.time,
                    self.evolution.time.last().copied().unwrap_or(f64::NAN)
                ),
            )
            .field("ψ-acc", &self.final_state.xacc)
            .field("θ-acc", &self.final_state.yacc)
            .field("Status", &self.status)
            .field("Steps taken", &self.evolution.steps_taken)
            .field("Steps stored", &self.evolution.time.len())
            .field("Initial energy", &self.initial_energy)
            .field("Final energy  ", &self.final_energy)
            .field("Duration", &self.evolution.duration)
            .finish()
    }
}
