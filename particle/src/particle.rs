use std::time::Duration;
use std::time::Instant;

use config::*;
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};

use crate::state::Display;
use crate::{Evolution, MappingParameters, PoincareSection, Solver, State};
use crate::{check_accuracy, map_integrate};

use crate::MagneticMoment;
use crate::{Distance, Flux, ParticleError, Radians, Result, Time};

use safe_unwrap::safe_unwrap;

/// A set of a Particle's intial conditions.
#[derive(Clone, Debug)]
pub struct InitialConditions {
    /// The initial time.
    pub time0: Time,
    /// The initial `θ` angle.
    pub theta0: Radians,
    /// The intial poloidal magnetic flux `ψp`.
    pub psip0: Flux,
    /// The initial parallel gyro radius `ρ`.
    pub rho0: Distance,
    /// The `ζ` angle.
    pub zeta0: Flux,
    /// The magnetic moment `μ`.
    pub mu: MagneticMoment,
}

#[derive(Debug, Clone, Default)]
pub enum IntegrationStatus {
    #[default]
    Initialized,
    Integrated,
    Escaped,
    TimedOut(Duration),
    InvalidIntersections,
    Failed {
        reason: String,
    },
}

/// Representation of a particle.
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
        let initial_state = State::from_initial(initial);
        let mut evolution = Evolution::with_capacity(EVOLUTION_INIT_CAPACITY);
        evolution.push_state(&initial_state);

        Self {
            initial_state,
            final_state: State::default(),
            status: IntegrationStatus::default(),
            evolution,
        }
    }

    /// Integrates the particle, storing the calculated orbit in [`Evolution`].
    pub fn integrate(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        currents: &Currents,
        perturbation: &Perturbation,
        t_eval: (Time, Time),
    ) -> Result<()> {
        self.initial_state
            .evaluate(qfactor, currents, bfield, perturbation)?;

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
                let result = state.evaluate(qfactor, currents, bfield, perturbation);
                self.evolution.push_state(&state);
                self.evolution.steps += 1;
                result
            };
            match res {
                Err(ParticleError::EqError(..)) => {
                    self.status = IntegrationStatus::Escaped;
                    break;
                }
                Err(err) => {
                    self.status = IntegrationStatus::Failed {
                        reason: format!("{:?}", err),
                    };
                    break;
                }
                Ok(_) => (),
            }

            if self.evolution.steps_taken() >= MAX_STEPS {
                self.status = IntegrationStatus::TimedOut(start.elapsed());
                break;
            }

            // Perform a step
            let mut solver = Solver::default();
            solver.init(&state);
            if let Err(err) = solver.start(dt, qfactor, bfield, currents, perturbation) {
                // This could only fail due to the solver's internal states' evaluate() calls.
                // However, this will be already caught at the start of the loop, even if the
                // initial state was invalid.
                unreachable!("{err}");
            };
            dt = solver.calculate_optimal_step(dt);
            state = solver.next_state(dt);
        }

        self.evolution.duration = start.elapsed();
        self.evolution.finish();
        self.final_state = state.into_evaluated(qfactor, currents, bfield, perturbation)?;
        Ok(())
    }

    /// Integrates the particle, storing its intersections with the Poincare surface defined by
    /// [`MappingParameters`].
    pub fn map(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        currents: &Currents,
        perturbation: &Perturbation,
        params: &MappingParameters,
    ) -> Result<()> {
        self.initial_state
            .evaluate(qfactor, currents, bfield, perturbation)?;
        self.status = IntegrationStatus::Integrated; // Will be overwritten in case of failure.
        let start = Instant::now();

        match map_integrate(self, qfactor, bfield, currents, perturbation, params) {
            Err(ParticleError::EqError(..)) => {
                self.status = IntegrationStatus::Escaped;
            }
            Err(ParticleError::TimedOut(..)) => {
                self.status = IntegrationStatus::TimedOut(start.elapsed());
            }
            Err(err) => {
                self.status = IntegrationStatus::Failed {
                    reason: format!("{:?}", err),
                };
            }
            Ok(_) => (),
        }

        let intersections = match params.section {
            PoincareSection::ConstZeta => &self.evolution.zeta,
            PoincareSection::ConstTheta => &self.evolution.theta,
        };
        if check_accuracy(intersections, MAP_THRESHOLD).is_err() {
            self.status = IntegrationStatus::InvalidIntersections;
        };

        self.evolution.duration = start.elapsed();
        self.evolution.finish();
        self.final_state = State {
            mu: self.initial_state.mu,
            time: safe_unwrap!("vec is non-empty", self.evolution.time.last().copied()),
            theta: safe_unwrap!("vec is non-empty", self.evolution.theta.last().copied()),
            psip: safe_unwrap!("vec is non-empty", self.evolution.psip.last().copied()),
            rho: safe_unwrap!("vec is non-empty", self.evolution.rho.last().copied()),
            zeta: safe_unwrap!("vec is non-empty", self.evolution.zeta.last().copied()),
            ..Default::default()
        }
        .into_evaluated(qfactor, currents, bfield, perturbation)?;

        Ok(())
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("ψ-acc", &self.final_state.xacc)
            .field("θ-acc", &self.final_state.yacc)
            .field("hcache", &self.final_state.hcache.first().or(None))
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
