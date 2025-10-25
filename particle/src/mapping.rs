use crate::Particle;
use crate::Point;
use crate::Result;
use crate::Solver;
use crate::State;

use equilibrium::{Bfield, Current, Perturbation, Qfactor};
use std::f64::consts::TAU;

/// Defines the surface of the Poincare section.
#[derive(Clone, Copy)]
pub enum PoincareSection {
    ConstTheta,
    ConstZeta,
}

/// Defines all the necessary parameters of a Poincare Map.
#[non_exhaustive]
#[derive(Clone, Copy)]
pub struct Mapping {
    /// The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    pub section: PoincareSection,
    /// The constant that defines the surface of section.
    pub alpha: f64,
    /// The number of interections to calculate.
    pub intersections: usize,
}

impl Mapping {
    pub fn new(section: PoincareSection, alpha: f64, intersections: usize) -> Self {
        // mod `alpha` to avoid modding it in every step
        Self {
            section,
            alpha: alpha % TAU,
            intersections,
        }
    }
}

/// Calculates the ζ=const intersections.
pub fn map_integrate(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    mapping: &Mapping,
) -> Result<()> {
    // Last two states of the particle.
    let mut state1 = particle.initial_state.clone(); // Already evaluated
    let mut state2: State;

    let mut dt = particle.config.rkf45_first_step;

    while particle.evolution.time.len() <= mapping.intersections {
        // Perform a step on the normal system.
        let mut solver = Solver::default();
        solver.init(&state1);
        solver.start(dt, qfactor, bfield, current, per)?;
        dt = solver.calculate_optimal_step(dt);
        state2 = solver.next_state(dt);
        state2.evaluate(qfactor, current, bfield, per)?;

        // Hénon's trick.
        // Depending on the PoincareSection, the independent variable becomes either `zeta` or
        // `theta`. Checking its value in every function and every loop has negligible performance
        // impact and produces much more readable code, instead of rewritting the same function
        // twice.
        if intersected(state1.zeta, state2.zeta, mapping.alpha) {
            let mod_state1 = calculate_mod_state1(&state1, &mapping.section);
            let dtau = calculate_mod_step(&state1, &state2, &mapping);
            let mod_state2 = calculate_mod_state2(qfactor, bfield, current, per, mod_state1, dtau)?;
            let intersection_state =
                calculate_intersection_state(qfactor, bfield, current, per, mapping, mod_state2)?;

            // Store the intersection state.
            particle
                .evolution
                .push_point(&Point::from_state(&intersection_state));

            // NOTE: Even after landing on the intersection, we must continue the integration from
            // state2. If we continue from the intersection state, a wrong sign change detection
            // will most likely occur, which causes the particle to get stuck. However, starting
            // from state2 is not a problem, since state2 was calculated from the solver, so it is
            // a valid state with a valid step size within the solver's tolerance.
        }
        // In both cases, continue from the next state.
        state1 = state2
    }
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
                ..*state1
            }
        }
    }
}

/// Calculates the step size dτ that brings mod_state1 on the intersection surface.
fn calculate_mod_step(state1: &State, state2: &State, mapping: &Mapping) -> f64 {
    // TODO: find a way to move the pole when the intersection angle is 0.
    match mapping.section {
        PoincareSection::ConstTheta => {
            let direction = (state2.theta - state1.theta).signum();
            direction * (mapping.alpha - state1.theta % TAU)
        }
        PoincareSection::ConstZeta => {
            let direction = (state2.zeta - state1.zeta).signum();
            direction * (mapping.alpha - state1.zeta % TAU)
        }
    }
}

/// Performs 1 step on the modified system (6) to calculate mod_state2, which sits exactly on the
/// intersection suface, **but corresponds to the modified system**.
fn calculate_mod_state2(
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    mod_state1: State,
    dtau: f64,
) -> Result<State> {
    let mut mod_solver = Solver::default();
    mod_solver.init(&mod_state1);
    mod_solver.start(dtau, qfactor, bfield, current, per)?;
    let mod_state2 = mod_solver.next_state(dtau);
    Ok(mod_state2)
}

/// Calculates the state of the original system exactly on the intersection surface, by converting
/// mod_state2 back on the original system.
fn calculate_intersection_state(
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    mapping: &Mapping,
    mod_state2: State,
) -> Result<State> {
    match mapping.section {
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
    .into_evaluated(qfactor, current, bfield, per)
}

// ====================================Common Functions===========================================

/// Checks when an angle has intersected with the surface at `angle`.
/// (source: seems to work)
fn intersected(old_angle: f64, new_angle: f64, surface_angle: f64) -> bool {
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() < 0.0
}

/// Checks if all the value diffs in the array are within the threshold.
pub fn check_accuracy(array: &[f64], threshold: f64) -> Result<()> {
    // array.iter().skip(1).for_each(|v| {
    //     dbg!(v % TAU);
    // });
    match array
        .windows(2)
        .skip(1) // Skip the starting point, since it is usually not on the intersection
        .all(|v| (v[1] % TAU - v[0] % TAU).abs() < threshold)
    {
        true => Ok(()),
        false => Err(crate::ParticleError::IntersectionError),
    }
}
