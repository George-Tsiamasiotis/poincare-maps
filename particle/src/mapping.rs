use crate::Particle;
use crate::ParticleError;
use crate::Radians;
use crate::Result;
use crate::Solver;
use crate::State;
use config::*;

use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use std::f64::consts::TAU;
use std::time::Duration;

/// Defines the surface of the Poincare section.
#[derive(Debug, Clone, Copy)]
pub enum PoincareSection {
    /// Defines a surface of xᵢ= θ.
    ConstTheta,
    /// Defines a surface of xᵢ= ζ.
    ConstZeta,
}

/// Defines all the necessary parameters of a Poincare Map.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct MappingParameters {
    /// The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    pub section: PoincareSection,
    /// The constant that defines the surface of section.
    pub alpha: Radians,
    /// The number of interections to calculate.
    pub intersections: usize,
}

impl MappingParameters {
    /// Creates a new [`MappingParameters`].
    pub fn new(section: PoincareSection, alpha: Radians, intersections: usize) -> Self {
        // mod `alpha` to avoid modding it in every step
        Self {
            section,
            alpha: alpha % TAU,
            intersections,
        }
    }
}

/// Calculates the PoincareSection=const intersections.
pub(crate) fn map_integrate(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    currents: &Currents,
    perturbation: &Perturbation,
    params: &MappingParameters,
) -> Result<()> {
    // Last two states of the particle.
    let mut state1 = particle.initial_state.clone(); // Already evaluated
    let mut state2: State;

    let mut dt = RKF45_FIRST_STEP;

    while particle.evolution.steps_stored() <= params.intersections {
        // Perform a step on the normal system.
        let mut solver = Solver::default();
        solver.init(&state1);
        solver.start(dt, qfactor, bfield, currents, perturbation)?;
        dt = solver.calculate_optimal_step(dt);
        state2 = solver.next_state(dt);
        state2.evaluate(qfactor, currents, bfield, perturbation)?;

        if particle.evolution.steps_taken() >= MAX_STEPS {
            return Err(ParticleError::TimedOut(Duration::default()));
        }

        // Hénon's trick.
        // Depending on the PoincareSection, the independent variable becomes either `zeta` or
        // `theta`. Checking its value in every function and every loop has negligible performance
        // impact and produces much more readable code, instead of rewritting the same function
        // twice.
        let (old_angle, new_angle) = match params.section {
            PoincareSection::ConstTheta => (state1.theta, state2.theta),
            PoincareSection::ConstZeta => (state1.zeta, state2.zeta),
        };
        if intersected(old_angle, new_angle, params.alpha) {
            let mod_state1 = calculate_mod_state1(&state1, &params.section);
            let dtau = calculate_mod_step(&state1, &state2, params);
            let mod_state2 =
                calculate_mod_state2(qfactor, bfield, currents, perturbation, mod_state1, dtau)?;
            let intersection_state = calculate_intersection_state(
                qfactor,
                bfield,
                currents,
                perturbation,
                params,
                mod_state2,
            )?;

            // Store the intersection state.
            particle.evolution.push_state(&intersection_state);

            // NOTE: Even after landing on the intersection, we must continue the integration from
            // state2. If we continue from the intersection state, a wrong sign change detection
            // will most likely occur, which causes the particle to get stuck. However, starting
            // from state2 is not a problem, since state2 was calculated from the solver, so it is
            // a valid state with a valid step size within the solver's tolerance.
        }
        // In both cases, continue from the next state.
        state1 = state2;
        particle.evolution.steps += 1;
    }
    particle.final_state = state1.into_evaluated(qfactor, currents, bfield, perturbation)?;
    Ok(())
}

/// Creates a [`State`] on the modified system (6) from `state1`, which is a state of the normal
/// system. The modified state corresponds to the first state on the modified system.
fn calculate_mod_state1(state1: &State, section: &PoincareSection) -> State {
    // Do not evaluate the state!
    match section {
        PoincareSection::ConstTheta => {
            let kappa = 1.0 / state1.theta_dot;
            let dt_dtheta = kappa;
            let dpsip_dtheta = kappa * state1.psip_dot;
            let drho_dtheta = kappa * state1.rho_dot;
            let dzeta_dtheta = kappa * state1.zeta_dot;
            State {
                time: state1.theta,
                theta: state1.time,
                theta_dot: dt_dtheta,
                psip_dot: dpsip_dtheta,
                rho_dot: drho_dtheta,
                zeta_dot: dzeta_dtheta,
                hcache: state1.hcache.clone(),
                ..*state1
            }
        }
        PoincareSection::ConstZeta => {
            let kappa = 1.0 / state1.zeta_dot;
            let dtheta_dzeta = kappa * state1.theta_dot;
            let dpsip_dzeta = kappa * state1.psip_dot;
            let drho_dzeta = kappa * state1.rho_dot;
            let dt_dzeta = kappa;
            State {
                time: state1.zeta,
                zeta: state1.time,
                theta_dot: dtheta_dzeta,
                psip_dot: dpsip_dzeta,
                rho_dot: drho_dzeta,
                zeta_dot: dt_dzeta,
                hcache: state1.hcache.clone(),
                ..*state1
            }
        }
    }
}

/// Calculates the step size dτ that brings mod_state1 on the intersection surface.
fn calculate_mod_step(state1: &State, state2: &State, params: &MappingParameters) -> f64 {
    // TODO: find a way to move the pole when the intersection angle is 0.
    match params.section {
        PoincareSection::ConstTheta => {
            let direction = (state2.theta - state1.theta).signum();
            direction * (params.alpha - state1.theta % TAU)
        }
        PoincareSection::ConstZeta => {
            let direction = (state2.zeta - state1.zeta).signum();
            direction * (params.alpha - state1.zeta % TAU)
        }
    }
}

/// Performs 1 step on the modified system (6) to calculate mod_state2, which sits exactly on the
/// intersection suface, **but corresponds to the modified system**.
fn calculate_mod_state2(
    qfactor: &Qfactor,
    bfield: &Bfield,
    currents: &Currents,
    perturbation: &Perturbation,
    mod_state1: State,
    dtau: f64,
) -> Result<State> {
    let mut mod_solver = Solver::default();
    mod_solver.init(&mod_state1);
    mod_solver.start(dtau, qfactor, bfield, currents, perturbation)?;
    let mod_state2 = mod_solver.next_state(dtau);
    Ok(mod_state2)
}

/// Calculates the state of the original system exactly on the intersection surface, by converting
/// mod_state2 back on the original system.
fn calculate_intersection_state(
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Currents,
    perturbation: &Perturbation,
    params: &MappingParameters,
    mod_state2: State,
) -> Result<State> {
    match params.section {
        PoincareSection::ConstTheta => {
            let kappa = 1.0;
            let dt_dt = kappa;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dzeta_dt = kappa * mod_state2.zeta_dot;
            State {
                time: mod_state2.theta,
                theta: mod_state2.time,
                theta_dot: dt_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dzeta_dt,
                mu: mod_state2.mu,
                ..mod_state2
            }
        }
        PoincareSection::ConstZeta => {
            let kappa = 1.0;
            let dtheta_dt = kappa * mod_state2.theta_dot;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dt_dt = kappa;
            State {
                time: mod_state2.zeta,
                zeta: mod_state2.time,
                theta_dot: dtheta_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dt_dt,
                mu: mod_state2.mu,
                ..mod_state2
            }
        }
    }
    .into_evaluated(qfactor, current, bfield, perturbation)
}

// ====================================Common Functions===========================================

/// Checks when an angle has intersected with the surface at `angle`.
/// (source: seems to work)
fn intersected(old_angle: Radians, new_angle: Radians, surface_angle: Radians) -> bool {
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    // NOTE: Use `<=` here for the case `surface_angle == 0`, since the sine of angles very close
    // to 0 (but not very close to 2π, 4π, ...)  return exactly 0.0.
    ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() <= 0.0
}

/// Checks if all the value diffs in the array are within the threshold.
pub(crate) fn check_accuracy(array: &[Radians], threshold: Radians) -> Result<()> {
    // array.iter().skip(1).for_each(|v| {
    //     dbg!(v % TAU);
    // });
    match array
        .windows(2)
        .skip(1) // Skip the starting point, since it is usually not on the intersection
        .all(|v| (v[1] - v[0]).abs() - TAU < threshold)
    {
        true => Ok(()),
        false => Err(crate::ParticleError::IntersectionError),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use std::f64::consts::{PI, TAU};

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
        // VERY corner case since all arguements have a sine of 0.0. Not worth working around. If
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

        assert!(check_accuracy(&ok1, MAP_THRESHOLD).is_ok());
        assert!(check_accuracy(&ok2, MAP_THRESHOLD).is_ok());
        assert!(check_accuracy(&ok3, MAP_THRESHOLD).is_ok());

        let not_ok1 = [
            0.0 * TAU,
            1.0 * TAU,
            // 2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let not_ok2 = [
            100.0,
            0.0 * TAU,
            1.0 * TAU,
            // 2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let not_ok3 = [
            100.0,
            1.0 * TAU,
            2.0 * TAU,
            3.0 * TAU + 1.0,
            4.0 * TAU,
            5.0 * TAU,
        ];

        assert!(check_accuracy(&not_ok1, MAP_THRESHOLD).is_err());
        assert!(check_accuracy(&not_ok2, MAP_THRESHOLD).is_err());
        assert!(check_accuracy(&not_ok3, MAP_THRESHOLD).is_err());
    }
}
