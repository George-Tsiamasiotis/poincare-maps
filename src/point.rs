use crate::State;

/// A point in the configuration space.
#[derive(Default, Debug, Clone)]
pub struct Point {
    pub time: f64,
    pub theta: f64,
    pub psip: f64,
    pub rho: f64,
    pub zeta: f64,
    pub mu: f64,
    pub ptheta: f64,
    pub pzeta: f64,
    pub psi: f64,
}

impl Point {
    pub fn to_state(&self) -> State {
        State {
            time: self.time,
            theta: self.theta,
            psip: self.psip,
            rho: self.rho,
            zeta: self.zeta,
            mu: self.mu,
            ..Default::default()
        }
    }
}
