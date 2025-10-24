use std::time::Instant;

use crate::Evolution;
use crate::ParticleError;
use crate::Point;
use crate::Result;
use crate::Solver;
use crate::State;
use crate::config::Config;
use crate::get_config;
use crate::state::Display;
use equilibrium::{Bfield, Current, Perturbation, Qfactor};

#[derive(Debug, Clone, Default)]
pub enum IntegrationStatus {
    #[default]
    Initialized,
    Integrated,
    Escaped,
    TimedOut,
    Failed,
}

#[derive(Clone)]
pub struct Particle {
    /// The initial [`State`] of the particle.
    pub initial_state: State,
    /// The final [`State`] of the particle.
    pub final_state: State,
    /// The [`Evolution`] time series of the particle.
    pub evolution: Evolution,
    /// Status of the particle's integration.
    pub status: IntegrationStatus,
    /// The global configuration parameters.
    config: Config,
}

impl Particle {
    /// Creates a new [`Particle`] from the initial conditions.
    pub fn new(t0: f64, theta0: f64, psip0: f64, rho0: f64, zeta0: f64, mu: f64) -> Self {
        let config = get_config();
        let point = Point::new(t0, theta0, psip0, rho0, zeta0, mu);
        let mut evolution = Evolution::with_capacity(config.evolution_init_capacity);
        evolution.push_point(&point);

        Self {
            initial_state: State::from_point(&point),
            final_state: State::default(),
            status: IntegrationStatus::default(),
            evolution,
            config,
        }
    }

    /// Integrates the particle, storing the calculated obrit in [`Evolution`].
    pub fn integrate(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        t_eval: (f64, f64),
    ) -> Result<()> {
        self.initial_state.evaluate(qfactor, current, bfield, per)?;

        // Tracks the state of the particle in each step. Also keeps the Accelerators' states.
        let mut state = self.initial_state.clone();
        self.status = IntegrationStatus::Integrated; // Will be overwritten in case of failure.
        let mut dt = self.config.rkf45_first_step;

        let start = Instant::now();
        while state.time <= t_eval.1 {
            // We still want to keep the particle, and also store its final state and final point,
            // even if the integration isn't properly completed.
            let res = {
                // Store the most recent state's point, including the intial and final points, even
                // if they are invalid.
                let result = state.evaluate(qfactor, current, bfield, per);
                self.evolution.push_point(&Point::from_state(&state));
                result
            };
            match res {
                Err(err) if matches!(&err, ParticleError::DomainError(..)) => {
                    self.status = IntegrationStatus::Escaped;
                    break;
                }
                Err(err) => {
                    // Should be unreachable..
                    eprintln!("{err}");
                    self.status = IntegrationStatus::Failed;
                    break;
                }
                Ok(_) => (),
            }

            if self.evolution.time.len() == self.config.max_steps {
                self.status = IntegrationStatus::TimedOut;
                return Err(ParticleError::TimedOut(start.elapsed()));
            }

            // Perform a step
            let mut solver = Solver::default();
            solver.init(&state);
            match solver.start(dt, qfactor, bfield, current, per) {
                Ok(_) => (),
                Err(err) => {
                    // This could only fail due to the solver's internal states' evaluate() calls.
                    // However, this will be already caught at the start of the loop, even if the
                    // initial state was invalid.
                    eprintln!("{err}");
                    unreachable!()
                }
            };
            dt = solver.calculate_optimal_step(dt);
            state = solver.next_state(dt);
        }

        self.evolution.duration = start.elapsed();
        self.final_state = state.into_evaluated(qfactor, current, bfield, per)?;
        Ok(())
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("ψ-acc", &self.final_state.xacc)
            .field("θ-acc", &self.final_state.yacc)
            .field("Initial", &Display::from_state(&self.initial_state))
            .field(
                "Initial parallel energy",
                &self.initial_state.parallel_energy(),
            )
            .field(
                "Intial perpendicular energy",
                &self.initial_state.perpendicular_energy(),
            )
            .field("Initial energy", &self.initial_state.energy())
            .field("Final energy  ", &self.final_state.energy())
            .field("Status", &self.status)
            .field("Evolution", &self.evolution)
            .finish()
    }
}
