#[cfg_attr(not(feature = "rkf45"), path = "rk45.rs")]
#[cfg_attr(feature = "rkf45", path = "rkf45.rs")]
mod rk;

use crate::State;
pub(crate) use rk::Solver;

/// Common to both solvers.
fn calculate_k1(solver: &mut Solver) {
    solver.k1 = [
        solver.state1.theta_dot,
        solver.state1.psip_dot,
        solver.state1.rho_dot,
        solver.state1.zeta_dot,
    ]
}

/// Common to both solvers.
fn next_state(solver: &mut Solver, h: f64) -> State {
    solver.next.t = solver.state1.t + h;
    solver.next.theta = solver.state1.theta + h * solver.weights[0];
    solver.next.psip = solver.state1.psip + h * solver.weights[1];
    solver.next.rho = solver.state1.rho + h * solver.weights[2];
    solver.next.zeta = solver.state1.zeta + h * solver.weights[3];
    solver.next.mu = solver.state1.mu;
    solver.next.clone()
}
