use std::f64::consts::TAU;

use crate::Particle;
use crate::Point;
use crate::Result;
use crate::Solver;
use crate::State;

use equilibrium::{Bfield, Current, Perturbation, Qfactor};

/// Defines the surface of the Poincare section.
pub enum PoincareSection {
    ConstTheta,
    ConstZeta,
}

/// Defines all the necessary parameters of a Poincare Map.
pub struct Mapping {
    /// The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    pub section: PoincareSection,
    /// The constant that defines the surface of section.
    pub alpha: f64,
    /// The number of interections to calculate.
    pub intersections: usize,
}

/// Calculates the ζ=const intersections.
pub fn zeta_map(
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
    let alpha_mod = mapping.alpha % TAU;

    while particle.evolution.time.len() <= mapping.intersections {
        // Perform a step on the normal system.
        let mut solver = Solver::default();
        solver.init(&state1);
        solver.start(dt, qfactor, bfield, current, per)?;
        dt = solver.calculate_optimal_step(dt);
        state2 = solver.next_state(dt);
        state2.evaluate(qfactor, current, bfield, per)?;

        // Hénon's trick
        if intersected(state1.zeta, state2.zeta, mapping.alpha) {
            // Switch to system (6).
            let kappa = 1.0 / state1.zeta_dot;
            let dtheta_dzeta = kappa * state1.theta_dot;
            let dpsip_dzeta = kappa * state1.psip_dot;
            let drho_dzeta = kappa * state1.rho_dot;
            let dt_dzeta = kappa;

            // First state in the modified system (6).
            let mod_state1 = State {
                time: state1.zeta,
                zeta: state1.time,
                theta_dot: dtheta_dzeta,
                psip_dot: dpsip_dzeta,
                rho_dot: drho_dzeta,
                zeta_dot: dt_dzeta,
                ..state1
            }; // Do not evaluate

            // Modified system step
            let direction = (state2.zeta - state1.zeta).signum();
            // TODO: find a way to move the pole when the intersection angle is 0.
            let dzeta = direction * (alpha_mod - state1.zeta % TAU);

            // Perform 1 step on the modified system (6). Do not use optimal step here.
            let mut mod_solver = Solver::default();
            mod_solver.init(&mod_state1);
            mod_solver.start(dzeta, qfactor, bfield, current, per)?;
            let mod_state2 = mod_solver.next_state(dzeta); // Do not evaluate

            // Switch back to the original system, which will give us the state at the surface.
            let kappa = 1.0;
            let dtheta_dt = kappa * mod_state2.theta_dot;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dt_dt = kappa;
            let mut intersection_state = State {
                time: mod_state2.zeta,
                zeta: mod_state2.time,
                theta_dot: dtheta_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dt_dt,
                mu: mod_state2.mu,
                ..mod_state2
            };
            intersection_state.evaluate(qfactor, current, bfield, per)?;

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

/// Calculates the θ=const intersections
#[allow(unused_variables)]
pub fn theta_map(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    mapping: &Mapping,
) -> Result<()> {
    todo!()
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
