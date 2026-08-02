//! Integration of a [`Particle`] for a specific time interval.

use std::time::Instant;

use dexter_equilibrium::Equilibrium;

use crate::particle::{IntegrationCaches, Particle};
use crate::solve::{SolverParams, Stepper};
use crate::state::GCState;

use super::IntegrationStatus;

// ===============================================================================================

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`IntegrationStatus`] variant for each possible error.
pub(super) fn integrate(
    particle: &mut Particle,
    equilibrium: &Equilibrium,
    teval: (f64, f64),
    solver_params: &SolverParams,
) {
    // =============== Setup

    let start = Instant::now();
    particle.evolution.reset();
    let mut caches = IntegrationCaches {
        mode_caches: equilibrium.perturbation.generate_caches(),
        ..Default::default()
    };

    // Return early if the initial flux happens to be exactly 0.0 or out of bounds.
    if particle.initial_conditions().flux0.value() == 0.0 {
        particle.integration_status = IntegrationStatus::OutOfBoundsInitialization;
        return;
    }
    if particle.initial_conditions.finalize(equilibrium).is_err() {
        particle.integration_status = IntegrationStatus::InvalidInitialConditions;
        return;
    }
    let Ok(mut state1) = GCState::new(&particle.initial_conditions, equilibrium, &mut caches)
    else {
        particle.integration_status = IntegrationStatus::OutOfBoundsInitialization;
        return;
    };

    particle.initial_energy = Some(state1.energy());
    let mut state2: GCState;
    let mut dt = solver_params.first_step;

    // =============== Main loop

    loop {
        if particle.evolution.steps_taken == solver_params.max_steps {
            particle.integration_status = IntegrationStatus::TimedOut(start.elapsed());
            break;
        }
        if particle.evolution.tf().is_some_and(|t| t > teval.1) {
            particle.integration_status = IntegrationStatus::Integrated;
            break;
        }

        // Perform a step
        let mut stepper = Stepper::new(&state1);
        state2 = if let Ok(state) = stepper
            .start(dt, equilibrium, &mut caches)
            .inspect(|_| dt = stepper.calculate_optimal_step(dt, solver_params))
            .and_then(|_| stepper.next_state(dt, equilibrium, &mut caches))
        {
            state
        } else {
            // `start()` and `next_state()` can only fail if an evaluation is out of bounds.
            particle.integration_status = IntegrationStatus::Escaped;
            break;
        };

        // Store and continue
        particle.evolution.push_state(&state1);
        particle.evolution.steps_taken += 1;
        state1 = state2;
    }

    // =============== Finalize

    particle.evolution.duration = start.elapsed();
    particle.final_energy = Some(state1.energy());
    particle.evolution.finish();
    particle.store_caches(caches);
}
