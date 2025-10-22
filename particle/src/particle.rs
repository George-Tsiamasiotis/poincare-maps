use std::time::{Duration, Instant};

use crate::Result;
use crate::get_config;
use crate::solver::Solver;
use crate::{InitialConditions, State};
use equilibrium::{Bfield, Current, Perturbation, Qfactor};

#[derive(Debug, Clone)]
pub enum ParticleStatus {
    Initialized(()),
    Integrated(()),
    Escaped(()),
    TimedOut(()),
}

#[derive(Clone)]
pub struct Particle {
    /// The initial (θ, ψ_p, ρ, ζ, μ) of the particle.
    pub initial: InitialConditions,
    /// The current state of the particle.
    pub state: State,
    /// The calculated evaluation times.
    pub t: Vec<f64>,
    /// The calculated θ values.
    pub theta: Vec<f64>,
    /// The calculated ψ_p values.
    pub psip: Vec<f64>,
    /// The calculated ρ values.
    pub rho: Vec<f64>,
    /// The calculated ζ values.
    pub zeta: Vec<f64>,
    /// The calculated Pζ values.
    pub pzeta: Vec<f64>,
    /// The calculated Pθ values.
    pub ptheta: Vec<f64>,
    /// The calculated ψ values.
    pub psi: Vec<f64>,
    /// The starting energy of the particle.
    pub initial_energy: f64,
    /// The final energy of the particle.
    pub final_energy: f64,
    /// The total compute time.
    pub calculation_time: Duration,
    /// The total number of integration steps.
    pub steps_taken: usize,
    /// Status about the particle's integration.
    pub status: ParticleStatus,
}

impl Particle {
    pub fn run_ode(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        t_eval: (f64, f64),
    ) -> Result<()> {
        self.state.evaluate(qfactor, current, bfield, per)?;
        self.initial_energy = self.state.energy();

        let mut h = get_config().rkf45_first_step;

        let start = Instant::now();

        while self.state.t < t_eval.1 {
            let mut solver = Solver::default();
            solver.init(&self.state);
            solver.start(h, qfactor, bfield, current, per)?;
            h = solver.calculate_optimal_step(h);
            self.state = solver.next_state(h);
            self.state.evaluate(qfactor, current, bfield, per)?;
            self.update_vecs();
            self.steps_taken += 1;
        }

        self.calculation_time = start.elapsed();
        self.status = ParticleStatus::Integrated(());
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
            .field("θ-acc", &self.state.yacc)
            .field("Status", &self.status)
            .field("Steps taken", &self.steps_taken)
            .field("Steps stored", &self.t.len())
            .field("Initial energy", &self.initial_energy)
            .field("Final energy  ", &self.final_energy)
            .field("Time", &self.calculation_time)
            .finish()
    }
}
