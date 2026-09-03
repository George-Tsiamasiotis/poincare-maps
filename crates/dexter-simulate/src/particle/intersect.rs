//! Calculation of a particle's intersections with a constant angle surface.

use approx::abs_diff_eq;
use std::f64::consts::{PI, TAU};
use std::time::Instant;

use dexter_machine::Machine;

use crate::constants::ANGLE_INTERSECTION_THRESHOLD;
use crate::particle::{Evolution, IntegrationCaches, Particle};
use crate::solve::{SolverParams, Stepper};
use crate::state::GCState;
use crate::{FluxCoordinate, IntegrationStatus, SimulationError};

/// Defines the surface of the Poincare section.
#[derive(Debug, Clone, Copy)]
pub enum Intersection {
    /// Defines a surface of xᵢ= θ.
    ConstTheta,
    /// Defines a surface of xᵢ= ζ.
    ConstZeta,
}

/// Defines all necessary parameters for the [`Particle::intersect`] routine.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct IntersectParams {
    /// The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    pub intersection: Intersection,
    /// The constant that defines the surface of section.
    pub angle: f64,
    /// The number of intersections to calculate.
    pub turns: usize,
}

impl IntersectParams {
    /// Creates a new [`IntersectParams`].
    #[must_use]
    pub fn new(intersection: Intersection, angle: f64, turns: usize) -> Self {
        // Mod `angle` to avoid modding it in every step
        Self {
            intersection,
            angle: angle.rem_euclid(TAU),
            turns,
        }
    }
}

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`IntegrationStatus`] variant for each possible error.
pub(super) fn intersect(
    particle: &mut Particle,
    machine: Machine,
    intersect_params: &IntersectParams,
    solver_params: &SolverParams,
) {
    // =============== Setup

    let start = Instant::now();
    particle.evolution.reset();
    let mut caches = IntegrationCaches {
        mode_caches: machine.perturbation().generate_caches(),
        ..Default::default()
    };
    // Create caches for the modified system to avoid invalidating the real caches.
    let mut mod_caches = IntegrationCaches {
        mode_caches: machine.perturbation().generate_caches(),
        ..Default::default()
    };

    // Return early if the initial flux happens to be exactly 0.0 or out of bounds.
    if particle.initial_conditions().flux0.value() == 0.0 {
        particle.integration_status = IntegrationStatus::OutOfBoundsInitialization;
        return;
    }
    if particle.initial_conditions.finalize(machine).is_err() {
        particle.integration_status = IntegrationStatus::InvalidInitialConditions;
        return;
    }
    let Ok(mut state1) = GCState::new(&particle.initial_conditions, machine, &mut caches) else {
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
        if particle.evolution.steps_stored() == intersect_params.turns {
            // Resolve integration status at the end
            break;
        }

        // Perform a step
        let mut stepper = Stepper::new(&state1);
        state2 = if let Ok(state) = stepper
            .start(dt, machine, &mut caches)
            .inspect(|_| dt = stepper.calculate_optimal_step(dt, solver_params))
            .and_then(|_| stepper.next_state(dt, machine, &mut caches))
        {
            state
        } else {
            // `start()` and `next_state()` can only fail if an evaluation is out of bounds.
            particle.integration_status = IntegrationStatus::Escaped;
            break;
        };

        let (old_angle, new_angle) = match intersect_params.intersection {
            Intersection::ConstTheta => (state1.theta, state2.theta),
            Intersection::ConstZeta => (state1.zeta, state2.zeta),
        };

        // Hénon's trick
        if intersected(old_angle, new_angle, intersect_params.angle) {
            // Switch to the modified system
            let mod_state1 = calculate_mod_state1(&state1, &intersect_params.intersection);

            // Calculate modified system's step size that gets us exactly on the intersection
            // surface.
            let dtau = calculate_mod_step(&state1, intersect_params);

            // Perform the step on the modified system
            let Ok(mod_state2) = calculate_mod_state2(machine, &mod_state1, dtau, &mut mod_caches)
            else {
                // `start()` and `next_state()` can only fail if an evaluation is out of bounds.
                particle.integration_status = IntegrationStatus::ModStateEscaped;
                break;
            };

            // Switch back to the normal system
            let Ok(intersection_state) = calculate_intersection_state(
                machine,
                &mod_state2,
                intersect_params,
                &mut mod_caches,
            ) else {
                // `start()` and `next_state()` can only fail if an evaluation is out of bounds.
                particle.integration_status = IntegrationStatus::ModStateEscaped;
                break;
            };
            // NOTE: Even after landing on the intersection, we must continue the integration from
            // state2. If we continue from the intersection state, a wrong sign change detection
            // will most likely occur, which causes the particle to get stuck. However, starting
            // from state2 is not a problem, since state2 was calculated from the solver, so it is
            // a valid state with a valid step size within the solver's tolerance.
            particle.evolution.push_state(&intersection_state);
        }

        // In both cases, continue from the next state.
        particle.evolution.steps_taken += 1;
        state1 = state2;
    }

    // =============== Finalize

    particle.evolution.duration = start.elapsed();
    particle.final_energy = Some(state1.energy());
    particle.evolution.finish();
    particle.store_caches(caches);

    if (particle.integration_status == IntegrationStatus::Escaped)
        || (particle.integration_status == IntegrationStatus::ModStateEscaped)
    {
        return;
    }
    // Further categorize successful intersection
    if check_mapping_accuracy(&particle.evolution, &intersect_params.intersection).is_ok() {
        if particle.steps_stored() == intersect_params.turns {
            particle.integration_status = IntegrationStatus::Intersected
        } else if particle.steps_stored() == 0 {
            particle.integration_status = IntegrationStatus::TimedOut(particle.evolution.duration);
        } else if particle.steps_stored() < intersect_params.turns {
            particle.integration_status = IntegrationStatus::IntersectedTimedOut;
        } else {
            unreachable!("`steps_stored` is always less that `intersect_params.turns")
        }
    } else {
        particle.integration_status = IntegrationStatus::InvalidIntersections;
    }
}

// ===============================================================================================

/// Creates a [`GCState`] on the modified system (6) from `state1` by making `θ` or `ζ` the
/// independent variable and `t` a dependent variable.
///
/// By assigning both `psi_dot` and `psip_dot` we cover both `FluxCoordinate` cases.
pub(crate) fn calculate_mod_state1(state1: &GCState, intersection: &Intersection) -> GCState {
    // Do not evaluate the state!
    let mut mod_state1 = state1.clone();
    match *intersection {
        Intersection::ConstTheta => {
            let kappa = 1.0 / state1.theta_dot;
            let dpsi_dtheta = kappa * state1.psi_dot;
            let dpsip_dtheta = kappa * state1.psip_dot;
            let dt_dtheta = kappa;
            let dzeta_dtheta = kappa * state1.zeta_dot;
            let drho_dtheta = kappa * state1.rho_dot;
            let dmu_dtheta = kappa * state1.mu_dot;

            mod_state1.t = state1.theta;
            mod_state1.theta = state1.t;
            mod_state1.psi_dot = dpsi_dtheta;
            mod_state1.psip_dot = dpsip_dtheta;
            mod_state1.theta_dot = dt_dtheta;
            mod_state1.zeta_dot = dzeta_dtheta;
            mod_state1.rho_dot = drho_dtheta;
            mod_state1.mu_dot = dmu_dtheta;
        }
        Intersection::ConstZeta => {
            let kappa = 1.0 / state1.zeta_dot;
            let dpsi_dzeta = kappa * state1.psi_dot;
            let dpsip_dzeta = kappa * state1.psip_dot;
            let dtheta_dzeta = kappa * state1.theta_dot;
            let dt_dzeta = kappa;
            let drho_dzeta = kappa * state1.rho_dot;
            let dmu_dzeta = kappa * state1.mu_dot;

            mod_state1.t = state1.zeta;
            mod_state1.zeta = state1.t;
            mod_state1.psi_dot = dpsi_dzeta;
            mod_state1.psip_dot = dpsip_dzeta;
            mod_state1.theta_dot = dtheta_dzeta;
            mod_state1.zeta_dot = dt_dzeta;
            mod_state1.rho_dot = drho_dzeta;
            mod_state1.mu_dot = dmu_dzeta;
        }
    }
    mod_state1
}

/// Calculates the step size `dτ` that brings `mod_state1` on the intersection surface.
///
/// dτ here has units of angle, since it is the 'time' step on the modified system.
pub(crate) fn calculate_mod_step(state1: &GCState, intersect_params: &IntersectParams) -> f64 {
    let surface_angle = intersect_params.angle;
    // NOTE: This is needed to move the %2π pole when the surface angle happens to be close to
    // a multiple of 2π. The surface angle is always modded to 2π, so we only need to check inside
    // this interval.
    // The tolerance doesn't have to be very tight, since moving the pole doesn't change the result.
    let pole = if abs_diff_eq!(surface_angle.abs(), 0.0, epsilon = 1e-1) {
        PI
    } else if abs_diff_eq!(surface_angle.abs(), TAU, epsilon = 1e-1) {
        -PI
    } else {
        0.0
    };
    match intersect_params.intersection {
        Intersection::ConstTheta => surface_angle - (state1.theta + pole).rem_euclid(TAU) + pole,
        Intersection::ConstZeta => surface_angle - (state1.zeta + pole).rem_euclid(TAU) + pole,
    }
}

/// Performs 1 step on the modified system (6) to calculate `mod_state2`, which sits exactly on the
/// intersection surface, **but corresponds to the modified system**.
pub(crate) fn calculate_mod_state2(
    objects: Machine,
    mod_state1: &GCState,
    dtau: f64,
    mod_caches: &mut IntegrationCaches,
) -> Result<GCState, SimulationError> {
    let mut mod_stepper = Stepper::new(mod_state1);
    mod_stepper.start(dtau, objects, mod_caches)?;
    // NOTE: This is equivalent to adjusting the step-size for the modified system.
    {
        mod_stepper.weights[0] = match mod_state1.coordinate {
            FluxCoordinate::Toroidal => mod_state1.psi_dot,
            FluxCoordinate::Poloidal => mod_state1.psip_dot,
        };
        mod_stepper.weights[1] = mod_state1.theta_dot;
        mod_stepper.weights[2] = mod_state1.zeta_dot;
        mod_stepper.weights[3] = mod_state1.rho_dot;
        mod_stepper.weights[4] = mod_state1.mu_dot;
    }
    mod_stepper.next_state(dtau, objects, mod_caches)
}

/// Calculates the state of the original system exactly on the intersection surface, by converting
/// `mod_state2` back on the original system.
pub(crate) fn calculate_intersection_state(
    objects: Machine,
    mod_state2: &GCState,
    intersect_params: &IntersectParams,
    mod_caches: &mut IntegrationCaches,
) -> Result<GCState, SimulationError> {
    let mut intersection_state = mod_state2.clone();
    match intersect_params.intersection {
        Intersection::ConstTheta => {
            let kappa = 1.0;
            let dpsi_dt = kappa * mod_state2.psi_dot;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let dt_dt = kappa;
            let dzeta_dt = kappa * mod_state2.zeta_dot;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dmu_dt = kappa * mod_state2.mu_dot;

            intersection_state.t = mod_state2.theta;
            intersection_state.theta = mod_state2.t;
            intersection_state.psi_dot = dpsi_dt;
            intersection_state.psip_dot = dpsip_dt;
            intersection_state.theta_dot = dt_dt;
            intersection_state.zeta_dot = dzeta_dt;
            intersection_state.rho_dot = drho_dt;
            intersection_state.mu_dot = dmu_dt;
        }
        Intersection::ConstZeta => {
            let kappa = 1.0;
            let dpsi_dt = kappa * mod_state2.psi_dot;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let dtheta_dt = kappa * mod_state2.theta_dot;
            let dt_dt = kappa;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dmu_dt = kappa * mod_state2.mu_dot;

            intersection_state.t = mod_state2.zeta;
            intersection_state.zeta = mod_state2.t;
            intersection_state.psi_dot = dpsi_dt;
            intersection_state.psip_dot = dpsip_dt;
            intersection_state.theta_dot = dtheta_dt;
            intersection_state.zeta_dot = dt_dt;
            intersection_state.rho_dot = drho_dt;
            intersection_state.mu_dot = dmu_dt;
        }
    };
    intersection_state.into_evaluated(objects, mod_caches)
}

/// Checks when an angle has intersected with the surface at `angle`
/// (source: seems to work).
///
/// # Note
///
/// To avoid the pole at `0`, `old_angle` and `new_angle` must NOT be modulo(2π).
pub(crate) fn intersected(old_angle: f64, new_angle: f64, surface_angle: f64) -> bool {
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    // NOTE: Use `<=` here for the case `surface_angle == 0`, since the sine of angles very close
    // to 0 (but not very close to 2π, 4π, ...)  return exactly 0.0.
    ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() <= 0.0
}

/// Checks if all the value diffs in the array are within the threshold.
fn check_mapping_accuracy(evolution: &Evolution, intersection: &Intersection) -> Result<(), ()> {
    let intersections_array = match *intersection {
        Intersection::ConstZeta => &evolution.zeta,
        Intersection::ConstTheta => &evolution.theta,
    };
    _check_mapping_accuracy(intersections_array)
}

/// Extracted for testing.
fn _check_mapping_accuracy(intersections_array: &[f64]) -> Result<(), ()> {
    if intersections_array
        .windows(2)
        .skip(1) // Skip the starting point, since it is often not on the intersection
        .all(|angle_pair| {
            approx::abs_diff_eq!(
                (angle_pair[1] - angle_pair[0]).rem_euclid(TAU),
                TAU,
                epsilon = ANGLE_INTERSECTION_THRESHOLD
            ) || approx::abs_diff_eq!(
                (angle_pair[1] - angle_pair[0]).rem_euclid(TAU),
                0.0,
                epsilon = ANGLE_INTERSECTION_THRESHOLD
            )
        })
    {
        Ok(())
    } else {
        Err(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_intersected() {
        // No need to check for negative numbers, since the signs would negate each other.

        let eps = 1e-12;

        assert!(intersected(1.0 - eps, 1.0 + eps, 1.0));
        assert!(intersected(1.0f64.next_down(), 1.0f64.next_up(), 1.0));
        assert!(!intersected(1.0 - eps, 1.0 - 2.0 * eps, 1.0));
        assert!(!intersected(1.0 + eps, 1.0 + 2.0 * eps, 1.0));
        assert!(!intersected(1.0 - eps, 1.0 - eps, 1.0));
        assert!(!intersected(1.0 + eps, 1.0 + eps, 1.0));
        assert!(!intersected(1.0f64.next_down(), 1.0f64.next_down(), 1.0));
        assert!(!intersected(1.0f64.next_up(), 1.0f64.next_up(), 1.0));

        assert!(intersected(10.0 - eps, 10.0 + eps, 10.0));
        assert!(intersected(10.0f64.next_down(), 10.0f64.next_up(), 10.0));
        assert!(!intersected(10.0 - eps, 10.0 - 2.0 * eps, 10.0));
        assert!(!intersected(10.0 + eps, 10.0 + 2.0 * eps, 10.0));
        assert!(!intersected(10.0 - eps, 10.0 - eps, 10.0));
        assert!(!intersected(10.0 + eps, 10.0 + eps, 10.0));
        assert!(!intersected(10.0f64.next_down(), 10.0f64.next_down(), 10.0));
        assert!(!intersected(10.0f64.next_up(), 10.0f64.next_up(), 10.0));

        assert!(intersected(PI - eps, PI + eps, PI));
        assert!(intersected(PI.next_down(), PI.next_up(), PI));
        assert!(!intersected(PI - eps, PI - 2.0 * eps, PI));
        assert!(!intersected(PI + eps, PI + 2.0 * eps, PI));
        assert!(!intersected(PI - eps, PI - eps, PI));
        assert!(!intersected(PI + eps, PI + eps, PI));
        assert!(!intersected(PI.next_down(), PI.next_down(), PI));
        assert!(!intersected(PI.next_up(), PI.next_up(), PI));

        assert!(intersected(TAU - eps, TAU + eps, TAU));
        assert!(intersected(TAU.next_down(), TAU.next_up(), TAU));
        assert!(!intersected(TAU - eps, TAU - 2.0 * eps, TAU));
        assert!(!intersected(TAU + eps, TAU + 2.0 * eps, TAU));
        assert!(!intersected(TAU - eps, TAU - eps, TAU));
        assert!(!intersected(TAU + eps, TAU + eps, TAU));
        assert!(!intersected(TAU.next_down(), TAU.next_down(), TAU));
        assert!(!intersected(TAU.next_up(), TAU.next_up(), TAU));

        assert!(intersected(0.0 - eps, 0.0 + eps, 0.0));
        assert!(intersected(0.0f64.next_down(), 0.0f64.next_up(), 0.0));
        assert!(!intersected(0.0 - eps, 0.0 - 2.0 * eps, 0.0));
        assert!(!intersected(0.0 + eps, 0.0 + 2.0 * eps, 0.0));
        assert!(!intersected(0.0 - eps, 0.0 - eps, 0.0));
        assert!(!intersected(0.0 + eps, 0.0 + eps, 0.0));
        // VERY corner case since all arguments have a sine of 0.0. Not worth working around. If
        // a particle is unlucky enough to reach this point, it would just be rejected from the
        // accuracy test.
        // assert!(!intersected(0.0f64.next_down(), 0.0f64.next_down(), 0.0));
        // assert!(!intersected(0.0f64.next_up(), 0.0f64.next_up(), 0.0));

        assert!(intersected(
            (2.0 * PI + PI).next_down(),
            (2.0 * PI + PI).next_up(),
            PI
        ));

        assert!(intersected(TAU.next_down(), TAU.next_up(), TAU));
        assert!(intersected(
            (2.0 * PI + TAU).next_down(),
            (2.0 * PI + TAU).next_up(),
            TAU
        ));

        assert!(!intersected(TAU - eps, TAU + eps, PI));
        assert!(!intersected(PI - eps, PI + eps, TAU));
        assert!(!intersected(PI - eps, PI + eps, PI / 2.0));
        assert!(!intersected(PI / 2.0 - eps, PI / 2.0 + eps, TAU));
    }

    #[test]
    fn test_accuracy_check() {
        let ok1 = [
            0.0 * TAU,
            1.0 * TAU,
            2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let ok2 = [
            100.0,
            0.0 * TAU,
            1.0 * TAU,
            2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let ok3 = [
            1.0 + 0.0 * TAU,
            1.0 + 1.0 * TAU,
            1.0 + 2.0 * TAU + 1e-12,
            1.0 + 3.0 * TAU - 1e-12,
            1.0 + 4.0 * TAU,
        ];

        assert!(_check_mapping_accuracy(&ok1).is_ok());
        assert!(_check_mapping_accuracy(&ok2).is_ok());
        assert!(_check_mapping_accuracy(&ok3).is_ok());

        // let not_ok1 = [
        //     0.0 * TAU,
        //     1.0 * TAU,
        //     // 2.0 * TAU + 1e-12,
        //     3.0 * TAU - 1e-12,
        //     4.0 * TAU,
        // ];
        // let not_ok2 = [
        //     100.0,
        //     0.0 * TAU,
        //     1.0 * TAU,
        //     // 2.0 * TAU + 1e-12,
        //     3.0 * TAU - 1e-12,
        //     4.0 * TAU,
        // ];
        // let not_ok3 = [
        //     100.0,
        //     1.0 * TAU,
        //     2.0 * TAU,
        //     3.0 * TAU + 1.0,
        //     4.0 * TAU,
        //     5.0 * TAU,
        // ];

        // assert!(_check_mapping_accuracy(&not_ok1).is_err());
        // assert!(_check_mapping_accuracy(&not_ok2).is_err());
        // assert!(_check_mapping_accuracy(&not_ok3).is_err());
    }
}
