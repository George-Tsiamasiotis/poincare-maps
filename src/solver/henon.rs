use std::time::Duration;
use std::time::Instant;

use crate::solver::Solver;
use crate::Particle;
use crate::Result;
use crate::State;
use crate::{Bfield, Current, Qfactor};
use std::f64::consts::TAU;

use crate::solver::RKF45_FIRST_STEP;

pub(crate) fn run_henon(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    angle: &str,
    intersection: f64,
    turns: usize,
) -> Result<()> {
    use std::f64::consts::TAU;
    particle.state.evaluate(qfactor, current, bfield)?;
    particle.initial_energy = particle.state.energy();

    particle.calculation_time = Duration::ZERO;
    let start = Instant::now();

    match angle {
        "zeta" => {
            henon_zeta_loop(particle, qfactor, bfield, current, intersection, turns)?;
            assert!(particle
                .zeta
                .windows(2)
                .map(|w| (w[1] - w[0]).abs() - TAU <= 1e-9)
                .all(|b| b));
        }
        "theta" => {
            henon_theta_loop(particle, qfactor, bfield, current, intersection, turns)?;
            assert!(particle
                .theta
                .windows(2)
                .map(|w| (w[1] - w[0]).abs() - TAU <= 1e-9)
                .all(|b| b));
        }
        _ => return Err(crate::MapError::InvalidAngle),
    };

    particle.calculation_time = start.elapsed();
    particle.shrink_vecs();

    particle.final_energy = particle.state.energy();

    Ok(())
}

pub(crate) fn henon_zeta_loop(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    intersection: f64,
    turns: usize,
) -> Result<()> {
    let mut h = RKF45_FIRST_STEP;

    while particle.zeta.len() <= turns {
        let (old_state, next_state) = get_step_stages(particle, qfactor, bfield, current, &mut h)?;

        if intersected(old_state.zeta, next_state.zeta, intersection) {
            // Setup new system (6) and (9)
            let kappa = 1.0 / old_state.zeta_dot;
            let dtheta_dzeta = kappa * old_state.theta_dot;
            let dpsip_dzeta = kappa * old_state.psip_dot;
            let drho_dzeta = kappa * old_state.rho_dot;
            let dt_dzeta = kappa;

            let old_mod_state = State {
                t: old_state.zeta,
                zeta: old_state.t,
                theta_dot: dtheta_dzeta,
                psip_dot: dpsip_dzeta,
                rho_dot: drho_dzeta,
                zeta_dot: dt_dzeta,
                ..old_state
            };

            let direction = (next_state.t - old_state.t).signum();
            // NOTE: unsure about where tho MOD should happen
            let dzeta = direction * (intersection - old_state.zeta % TAU);

            let next_mod_state =
                calculate_next_mod_state(qfactor, bfield, current, old_mod_state, dzeta)?;

            let kappa = 1.0;
            let dtheta_dt = kappa * next_mod_state.theta_dot;
            let dpsip_dt = kappa * next_mod_state.psip_dot;
            let drho_dt = kappa * next_mod_state.rho_dot;
            let dt_dt = kappa;
            let intersection_state = State {
                t: next_mod_state.zeta,
                zeta: next_mod_state.t,
                theta_dot: dtheta_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dt_dt,
                ..next_mod_state
            };

            store_intersection(
                particle,
                qfactor,
                bfield,
                current,
                &next_state,
                intersection_state,
            )?;
            h = RKF45_FIRST_STEP;
        } else {
            particle.state = next_state;
        }
    }
    Ok(())
}

pub(crate) fn henon_theta_loop(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    intersection: f64,
    turns: usize,
) -> Result<()> {
    let mut h = RKF45_FIRST_STEP;

    while particle.theta.len() <= turns {
        let (old_state, next_state) = get_step_stages(particle, qfactor, bfield, current, &mut h)?;

        if intersected(old_state.theta, next_state.theta, intersection) {
            // Setup new system (6) and (9)
            let kappa = 1.0 / old_state.theta_dot;
            let dt_dtheta = kappa;
            let dpsip_dtheta = kappa * old_state.psip_dot;
            let drho_dtheta = kappa * old_state.rho_dot;
            let dzeta_dtheta = kappa * old_state.zeta_dot;

            let old_mod_state = State {
                t: old_state.theta,
                theta: old_state.t,
                theta_dot: dt_dtheta,
                psip_dot: dpsip_dtheta,
                rho_dot: drho_dtheta,
                zeta_dot: dzeta_dtheta,
                ..old_state
            };

            let direction = (next_state.t - old_state.t).signum();
            // NOTE: unsure about where tho MOD should happen
            let dtheta = direction * (intersection - old_state.theta % TAU);

            let next_mod_state =
                calculate_next_mod_state(qfactor, bfield, current, old_mod_state, dtheta)?;

            let kappa = 1.0;
            let dt_dt = kappa;
            let dpsip_dt = kappa * next_mod_state.psip_dot;
            let drho_dt = kappa * next_mod_state.rho_dot;
            let dzeta_dt = kappa * next_mod_state.zeta_dot;
            let intersection_state = State {
                t: next_mod_state.theta,
                theta: next_mod_state.t,
                theta_dot: dt_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dzeta_dt,
                ..next_mod_state
            };

            store_intersection(
                particle,
                qfactor,
                bfield,
                current,
                &next_state,
                intersection_state,
            )?;
            h = RKF45_FIRST_STEP;
        } else {
            particle.state = next_state;
        }
    }
    Ok(())
}

// ====================================Common Functions===========================================

/// Returns the (normal) system's intial and next state, by performing 1 step from the initial
/// state.
fn get_step_stages(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    h: &mut f64,
) -> Result<(State, State)> {
    let mut solver = Solver::default();
    solver.init(&particle.state);
    solver.start(*h, qfactor, bfield, current)?;
    *h = solver.calculate_optimal_step(*h);
    let old_state = particle.state.clone();
    let mut next_state = solver.next_state(*h);
    next_state.evaluate(qfactor, current, bfield)?;
    Ok((old_state, next_state))
}

/// Checks when an angle has intersected with the surface at `angle`.
/// (source: seems to work)
fn intersected(old_angle: f64, new_angle: f64, surface_angle: f64) -> bool {
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() < 0.0
}

/// Performs 1 step on the modified system and returns the new modified state.
fn calculate_next_mod_state(
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    mod_state: State,
    dtau: f64,
) -> Result<State> {
    let mut solver = Solver::default();
    solver.init(&mod_state);
    solver.start(dtau, qfactor, bfield, current)?;
    let new_mod_state = solver.next_state(dtau);
    Ok(new_mod_state)
}

/// Evaluates and stores the State at the intersction surface in the particle.
fn store_intersection(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    next_state: &State,
    intersection_state: State,
) -> Result<()> {
    particle.state = intersection_state.clone();
    particle.state.evaluate(qfactor, current, bfield)?;
    particle.update_vecs();
    let temp_step = next_state.t - intersection_state.t;
    let mut solver = Solver::default();
    solver.init(&particle.state);
    solver.start(temp_step, qfactor, bfield, current)?;
    particle.state = solver.next_state(temp_step);
    particle.state.evaluate(qfactor, current, bfield)?;
    Ok(())
}
