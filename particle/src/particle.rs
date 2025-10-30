use std::time::Duration;
use std::time::Instant;

use crate::Evolution;
use crate::InitialConditions;
use crate::Mapping;
use crate::ParticleError;
use crate::PoincareSection;
use crate::Point;
use crate::Result;
use crate::Solver;
use crate::State;
use crate::check_accuracy;
use crate::map_integrate;
use crate::state::Display;
use config::*;
use equilibrium::{Bfield, Current, Perturbation, Qfactor};

#[derive(Debug, Clone, Default)]
pub enum IntegrationStatus {
    #[default]
    Initialized,
    Integrated,
    Escaped,
    TimedOut(Duration),
    InvalidIntersections,
    Failed {
        reason: Box<str>,
    },
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
}

impl Particle {
    /// Creates a new [`Particle`] from the initial conditions.
    pub fn new(initial: &InitialConditions) -> Self {
        let point = initial.to_point();
        let mut evolution = Evolution::with_capacity(EVOLUTION_INIT_CAPACITY);
        evolution.push_point(&point);

        Self {
            initial_state: State::from_point(&point),
            final_state: State::default(),
            status: IntegrationStatus::default(),
            evolution,
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
        let mut dt = RKF45_FIRST_STEP;

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
                Err(ParticleError::DomainError(..)) => {
                    self.status = IntegrationStatus::Escaped;
                    break;
                }
                Err(err) => {
                    self.status = IntegrationStatus::Failed { reason: err.into() };
                    break;
                }
                Ok(_) => (),
            }

            if self.evolution.time.len() == MAX_STEPS {
                self.status = IntegrationStatus::TimedOut(start.elapsed());
                break;
            }

            // Perform a step
            let mut solver = Solver::default();
            solver.init(&state);
            if let Err(err) = solver.start(dt, qfactor, bfield, current, per) {
                // This could only fail due to the solver's internal states' evaluate() calls.
                // However, this will be already caught at the start of the loop, even if the
                // initial state was invalid.
                unreachable!("{err}");
            };
            dt = solver.calculate_optimal_step(dt);
            state = solver.next_state(dt);
        }

        self.evolution.duration = start.elapsed();
        self.evolution.shrink_to_fit();
        self.final_state = state.into_evaluated(qfactor, current, bfield, per)?;
        Ok(())
    }

    /// Integrates the particle, storing its intersections with the Poincare surface defined by
    /// [`Mapping`].
    pub fn map(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        mapping: &Mapping,
    ) -> Result<()> {
        self.initial_state.evaluate(qfactor, current, bfield, per)?;
        self.status = IntegrationStatus::Integrated; // Will be overwritten in case of failure.
        let start = Instant::now();

        match map_integrate(self, qfactor, bfield, current, per, mapping) {
            Err(ParticleError::DomainError(..)) => {
                self.status = IntegrationStatus::Escaped;
            }
            Err(ParticleError::TimedOut(..)) => {
                self.status = IntegrationStatus::TimedOut(start.elapsed());
            }
            Err(err) => {
                self.status = IntegrationStatus::Failed { reason: err.into() };
            }
            Ok(_) => (),
        }

        let intersections = match mapping.section {
            PoincareSection::ConstZeta => &self.evolution.zeta,
            PoincareSection::ConstTheta => &self.evolution.theta,
        };
        if let Some(err) = check_accuracy(intersections, MAP_THRESHOLD).err() {
            self.status = IntegrationStatus::Failed { reason: err.into() };
        };

        self.evolution.duration = start.elapsed();
        self.evolution.shrink_to_fit();
        self.final_state = State {
            mu: self.initial_state.mu,
            time: self.evolution.time.last().copied().unwrap_or_default(),
            theta: self.evolution.theta.last().copied().unwrap_or_default(),
            psip: self.evolution.psip.last().copied().unwrap_or_default(),
            rho: self.evolution.rho.last().copied().unwrap_or_default(),
            zeta: self.evolution.zeta.last().copied().unwrap_or_default(),
            ..Default::default()
        }
        .into_evaluated(qfactor, current, bfield, per)?;

        Ok(())
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("ψ-acc", &self.final_state.xacc)
            .field("θ-acc", &self.final_state.yacc)
            .field("hcache", &self.final_state.hcache.first().unwrap())
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
