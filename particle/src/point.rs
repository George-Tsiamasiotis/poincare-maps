use crate::State;

/// Representation of a single pointing, containing the independent variables, but also the derived
/// toroidal flux and canonical momenta.
pub struct Point {
    pub time: f64,
    pub theta: f64,
    pub psip: f64,
    pub rho: f64,
    pub zeta: f64,
    pub mu: f64,
    pub psi: f64,
    pub ptheta: f64,
    pub pzeta: f64,
}

impl Point {
    /// Creates a new Point from the **independent** variables.
    pub fn new(time: f64, theta: f64, psip: f64, rho: f64, zeta: f64, mu: f64) -> Self {
        Self {
            time,
            theta,
            psip,
            rho,
            zeta,
            mu,
            ..Point::default()
        }
    }

    /// Creates a new Point from a state, copying both the independent and the derived variables.
    pub fn from_state(state: &State) -> Self {
        Self {
            time: state.time,
            theta: state.theta,
            psip: state.psip,
            rho: state.rho,
            zeta: state.zeta,
            mu: state.mu,
            psi: state.psi,
            ptheta: state.ptheta,
            pzeta: state.pzeta,
        }
    }
}

impl Default for Point {
    fn default() -> Self {
        Self {
            time: f64::NAN,
            theta: f64::NAN,
            psip: f64::NAN,
            rho: f64::NAN,
            zeta: f64::NAN,
            psi: f64::NAN,
            mu: f64::NAN,
            ptheta: f64::NAN,
            pzeta: f64::NAN,
        }
    }
}
