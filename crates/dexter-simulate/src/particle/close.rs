//! Integration of a [`Particle`] for a specific amound of full `θ-ψ` periods.

use std::f64::consts::TAU;
use std::time::Instant;

use approx::relative_eq;
use dexter_equilibrium::Equilibrium;

use crate::constants::{FLUX_REL_TOL, SHORT_CIRCUIT_FLUX_REL_TOL};
use crate::particle::intersect::{
    calculate_intersection_state, calculate_mod_state1, calculate_mod_state2, calculate_mod_step,
    intersected,
};
use crate::particle::{IntegrationCaches, Particle};
use crate::solve::{SolverParams, Stepper};
use crate::state::GCState;
use crate::{Frequencies, IntersectParams};

use super::{IntegrationStatus, Intersection};

// ===============================================================================================

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`IntegrationStatus`] variant for each possible error.
pub(super) fn close(
    particle: &mut Particle,
    equilibrium: &Equilibrium,
    periods: usize,
    solver_params: &SolverParams,
) {
    // =============== Setup

    let start = Instant::now();
    particle.evolution.reset();
    let mut caches = IntegrationCaches {
        mode_caches: equilibrium.perturbation.generate_caches(),
        ..Default::default()
    };
    let mut mod_caches = IntegrationCaches {
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
    let Ok(state0) = GCState::new(&particle.initial_conditions, equilibrium, &mut caches) else {
        particle.integration_status = IntegrationStatus::OutOfBoundsInitialization;
        return;
    };

    // `state0` has been evaluated on the initial point
    particle.initial_energy = Some(state0.energy());
    let mut state1 = state0.clone();

    let mut state2: GCState;
    let mut dt = solver_params.first_step;
    let mut closed_periods: usize = 0;

    // Phony intersection parameters to be used in the modified system. The `turns` field is
    // ignored.
    let intersect_params = &IntersectParams::new(Intersection::ConstTheta, state0.theta, 0);

    // =============== Main loop

    loop {
        if particle.evolution.steps_taken == solver_params.max_steps {
            particle.integration_status = match closed_periods {
                0 => IntegrationStatus::TimedOut(start.elapsed()),
                1.. => IntegrationStatus::ClosedPeriods(closed_periods),
            };
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

        // This may only fail if an evaluation of the modified system fails
        let Ok(close_period_check) = closed_period(
            equilibrium,
            &state0,
            &state1,
            &state2,
            intersect_params,
            &mut mod_caches,
        ) else {
            particle.integration_status = IntegrationStatus::ModStateEscaped;
            break;
        };

        // Avoid stopping the particle at the start of the integration. When the `stepping_method`
        // is `EnergyAdaptiveStep`, the step size grows almost exponentially, so 10 steps should be
        // enough for `ψ` to move sufficiently far from `ψ0`.
        if particle.steps_taken() > 10 && close_period_check {
            closed_periods += 1;
        };

        if closed_periods == periods {
            particle.integration_status = IntegrationStatus::ClosedPeriods(periods);
            particle.evolution.push_state(&state1);
            break;
        }

        // Store and continue
        particle.evolution.push_state(&state1);
        particle.evolution.steps_taken += 1;
        state1 = state2;
    }

    // =============== Hénon trick to close the period

    // When we are about to close the final period, we use Hénon's trick to take one more step
    // and land exactly on the initial point. This significantly improves the accuracy of the
    // calculations of the frequencies.
    if particle.integration_status == IntegrationStatus::ClosedPeriods(periods) {
        // Switch to the modified system
        let mod_state1 = calculate_mod_state1(&state1, &intersect_params.intersection);

        // Time step to accurately close the orbit
        let dtau = calculate_mod_step(&state1, intersect_params);

        // Perform the step on the modified system
        // Can only fail if an evaluation is out of bounds
        //
        // If `mod_state2` was calculated correctly, switch back to the normal system, calculate
        // `intersection_state` and store it
        match calculate_mod_state2(equilibrium, &mod_state1, dtau, &mut mod_caches) {
            Ok(mod_state2) => {
                match calculate_intersection_state(
                    equilibrium,
                    &mod_state2,
                    intersect_params,
                    &mut mod_caches,
                ) {
                    Ok(intersection_state) => {
                        particle.evolution.push_state(&intersection_state);
                        particle.integration_status =
                            IntegrationStatus::ClosedPeriods(closed_periods)
                    }
                    Err(_) => {
                        particle.integration_status = IntegrationStatus::ModStateEscaped;
                    }
                }
            }
            Err(_) => {
                particle.integration_status = IntegrationStatus::ModStateEscaped;
            }
        };
    }

    // =============== Finalize

    calculate_frequencies(particle);
    particle.evolution.duration = start.elapsed();
    particle.final_energy = Some(state1.energy());
    particle.evolution.finish();
    particle.store_caches(caches);
}

/// Checks if `state1` and `state2` stand 'left and right' of the period closing point. This
/// indicates that the particle reached its starting point.
///
/// The check is performed by taking a step on the modified system, which guarantees that `θ` will
/// be exactly `θ0`. Then, it checks if `ψ` is also close to `ψ0`, which indicates that the orbit is
/// closed.
///
/// # Errors
///
/// This method returns an `Err(())` if an evaluation on the modified system fails, which is the
/// only possible error.
fn closed_period(
    equilibrium: &Equilibrium,
    state0: &GCState,
    state1: &GCState,
    state2: &GCState,
    intersect_params: &IntersectParams,
    mod_caches: &mut IntegrationCaches,
) -> Result<bool, ()> {
    // PERF: Short-circuit the `ψ` check to avoid the more expensive checks
    if !relative_eq!(state0.psi, state1.psi, epsilon = SHORT_CIRCUIT_FLUX_REL_TOL) {
        return Ok(false);
    }
    if !intersected(state1.theta, state2.theta, state0.theta) {
        return Ok(false);
    }

    // Switch to the modified system
    let mod_state1 = calculate_mod_state1(state1, &intersect_params.intersection);
    // Time step to accurately close the orbit
    let dtau = calculate_mod_step(state1, intersect_params);

    // Perform the step on the modified system
    //
    // If `mod_state2` was calculated correctly, switch back to the normal system, calculate
    // `intersection_state` and check for closed period
    //
    // The two calculations can only fail if an evaluation of the modified system is out of bounds
    let Ok(mod_state2) = calculate_mod_state2(equilibrium, &mod_state1, dtau, mod_caches) else {
        return Err(());
    };
    let Ok(intersection_state) =
        calculate_intersection_state(equilibrium, &mod_state2, intersect_params, mod_caches)
    else {
        return Err(());
    };

    // Final `ψ-ψ0` and `θ` direction checks
    if !relative_eq!(state0.psi, intersection_state.psi, epsilon = FLUX_REL_TOL) {
        return Ok(false);
    }
    if state0.theta_dot.signum() != intersection_state.theta_dot.signum() {
        return Ok(false);
    }
    Ok(true)
}

/// Calculates the final `ωθ`, `ωζ` and `qkinetic` and stores them in the particle's
/// `frequencies` field.
///
/// Also takes account the number of periods closed.
fn calculate_frequencies(particle: &mut Particle) {
    let periodsf64 = match particle.integration_status {
        IntegrationStatus::ClosedPeriods(periods) => periods as f64,
        _ => return,
    };
    // Use `NaN`for particles with invalid initial conditions or out of bounds initialization.
    let t0 = particle.evolution.t.first().copied().unwrap_or(f64::NAN);
    let tf = particle.evolution.t.last().copied().unwrap_or(f64::NAN);
    let theta_period = tf - t0;

    let zeta0 = particle.evolution.zeta.first().copied().unwrap_or(f64::NAN);
    let zetaf = particle.evolution.zeta.last().copied().unwrap_or(f64::NAN);
    let dzeta = zetaf - zeta0;

    let omega_theta = TAU / (theta_period / periodsf64);
    let omega_zeta = dzeta / (theta_period / periodsf64);
    let qkinetic = omega_zeta / omega_theta;

    particle.frequencies = Frequencies {
        omega_theta: Some(omega_theta),
        omega_zeta: Some(omega_zeta),
        qkinetic: Some(qkinetic),
    }
}
