use std::time::Instant;

use crate::Evolution;
use crate::Particle;
use crate::PoincareParameters;
use crate::Result;
use crate::State;
use crate::Surface;
use crate::solver::Solver;
use crate::{Bfield, Current, Perturbation, Qfactor};
use std::f64::consts::TAU;

use crate::MAX_STEPS;
use crate::RKF45_FIRST_STEP;

pub(crate) fn run_henon(
    evolution: &mut Evolution,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    params: &PoincareParameters,
) -> Result<()> {
    let start = Instant::now();

    match params.surface {
        Surface::ConstZeta => {
            henon_zeta_loop(evolution, qfactor, bfield, current, per, params)?;
            check_accuracy(&evolution.zeta);
        }
        Surface::ConstTheta => {
            // henon_theta_loop(evolution, qfactor, bfield, current, per, params)?;
            // check_accuracy(&evolution.theta);
        }
    };

    evolution.duration = start.elapsed();
    evolution.shrink_to_fit();

    Ok(())
}

pub(crate) fn henon_zeta_loop(
    evolution: &mut Evolution,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    params: &PoincareParameters,
) -> Result<()> {
    let calculation_time = Instant::now();
    let mut h = RKF45_FIRST_STEP;

    while evolution.steps_stored < params.turns {
        let last_state = evolution.last_state();
        let next_state = get_next_stage(qfactor, bfield, current, &last_state, per, &mut h)?
            .into_evaluated(qfactor, current, bfield, per)?;

        evolution.steps_taken += 1;
        if evolution.steps_taken >= MAX_STEPS {
            return Err(crate::MapError::OrbitTimeout(
                calculation_time.elapsed(),
                MAX_STEPS,
            ));
        }

        if intersected(last_state.zeta, next_state.zeta, params.intersection) {
            // Setup new system (6) and (9)
            let kappa = 1.0 / last_state.zeta_dot;
            let dtheta_dzeta = kappa * last_state.theta_dot;
            let dpsip_dzeta = kappa * last_state.psip_dot;
            let drho_dzeta = kappa * last_state.rho_dot;
            let dt_dzeta = kappa;

            let old_mod_state = State {
                time: last_state.zeta,
                zeta: last_state.time,
                theta_dot: dtheta_dzeta,
                psip_dot: dpsip_dzeta,
                rho_dot: drho_dzeta,
                zeta_dot: dt_dzeta,
                ..last_state
            };

            let direction = (next_state.time - last_state.time).signum();
            // NOTE: unsure about where tho MOD should happen
            let mut dzeta = direction * (params.intersection - last_state.zeta % TAU);

            let next_mod_state = // Do not evaluate
                get_next_stage(qfactor, bfield, current, &old_mod_state, per, &mut dzeta)?;

            let kappa = 1.0;
            let dtheta_dt = kappa * next_mod_state.theta_dot;
            let dpsip_dt = kappa * next_mod_state.psip_dot;
            let drho_dt = kappa * next_mod_state.rho_dot;
            let dt_dt = kappa;
            let intersection_state = State {
                time: next_mod_state.zeta,
                zeta: next_mod_state.time,
                theta_dot: dtheta_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dt_dt,
                ..next_mod_state
            };

            // Store the Intersection state.
            evolution.push_point(
                intersection_state
                    .into_evaluated(qfactor, current, bfield, per)?
                    .as_point(),
            );
        } else {
            // Just keep going if no intersection was detected
            evolution.push_point(next_state.as_point());
        }
    }
    Ok(())
}

#[allow(dead_code)]
#[allow(unused_variables)]
pub(crate) fn henon_theta_loop(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    per: &Perturbation,
    intersection: f64,
    turns: usize,
) -> Result<()> {
    // let mut h = RKF45_FIRST_STEP;
    // let calculation_time = Instant::now();
    // while particle.evolution.theta.len() < turns {
    //     let (old_state, next_state) =
    //         get_next_stage(particle, qfactor, bfield, current, per, &mut h)?;
    //     particle.steps_taken += 1;
    //     if particle.steps_taken >= MAX_STEPS {
    //         return Err(crate::MapError::OrbitTimeout(
    //             calculation_time.elapsed(),
    //             MAX_STEPS,
    //         ));
    //     }
    //
    //     if intersected(old_state.theta, next_state.theta, intersection) {
    //         // Setup new system (6) and (9)
    //         let kappa = 1.0 / old_state.theta_dot;
    //         let dt_dtheta = kappa;
    //         let dpsip_dtheta = kappa * old_state.psip_dot;
    //         let drho_dtheta = kappa * old_state.rho_dot;
    //         let dzeta_dtheta = kappa * old_state.zeta_dot;
    //
    //         let old_mod_state = State {
    //             time: old_state.theta,
    //             theta: old_state.time,
    //             theta_dot: dt_dtheta,
    //             psip_dot: dpsip_dtheta,
    //             rho_dot: drho_dtheta,
    //             zeta_dot: dzeta_dtheta,
    //             ..old_state
    //         };
    //
    //         let direction = (next_state.time - old_state.time).signum();
    //         // NOTE: unsure about where tho MOD should happen
    //         let dtheta = direction * (intersection - old_state.theta % TAU);
    //
    //         let next_mod_state =
    //             calculate_next_mod_state(qfactor, bfield, current, per, old_mod_state, dtheta)?;
    //
    //         let kappa = 1.0;
    //         let dt_dt = kappa;
    //         let dpsip_dt = kappa * next_mod_state.psip_dot;
    //         let drho_dt = kappa * next_mod_state.rho_dot;
    //         let dzeta_dt = kappa * next_mod_state.zeta_dot;
    //         let intersection_state = State {
    //             time: next_mod_state.theta,
    //             theta: next_mod_state.time,
    //             theta_dot: dt_dt,
    //             psip_dot: dpsip_dt,
    //             rho_dot: drho_dt,
    //             zeta_dot: dzeta_dt,
    //             ..next_mod_state
    //         };
    //
    //         store_intersection(
    //             particle,
    //             qfactor,
    //             bfield,
    //             current,
    //             per,
    //             &next_state,
    //             intersection_state,
    //         )?;
    //         h = RKF45_FIRST_STEP;
    //     } else {
    //         particle.state = next_state;
    //     }
    // }
    Ok(())
}

// ====================================Common Functions===========================================

/// Returns the (normal) system's next state and the next optimal step, by performing 1 step from
/// the last state. The last state is **not** evaluated.
fn get_next_stage(
    qfactor: &Qfactor,
    bfield: &Bfield,
    current: &Current,
    last_state: &State,
    per: &Perturbation,
    h: &mut f64,
) -> Result<State> {
    let mut solver = Solver::default();
    solver.init(last_state);
    solver.start(*h, qfactor, bfield, current, per)?;
    *h = solver.calculate_optimal_step(*h);
    let next_state = solver.next_state(*h);
    Ok(next_state)
}

/// Checks when an angle has intersected with the surface at `angle`.
/// (source: seems to work)
fn intersected(old_angle: f64, new_angle: f64, surface_angle: f64) -> bool {
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() < 0.0
}

fn check_accuracy(angle_time_series: &[f64]) {
    assert!(
        angle_time_series
            .windows(2)
            .all(|w| (w[1] - w[0]).abs() - TAU <= 1e-9)
    );
}
