#[derive(Clone, Debug)]
pub struct InitialConditions {
    /// The initial time.
    pub t0: f64,
    /// The initial `θ` angle.
    pub theta0: f64,
    /// The intial poloidal magnetic flux `ψ_p`.
    pub psip0: f64,
    /// The initial parallel gyro radius `ρ`.
    pub rho0: f64,
    /// The `ζ` angle.
    pub zeta0: f64,
    /// The magnetic moment `μ`.
    pub mu: f64,
}

impl InitialConditions {
    pub fn new(t0: f64, theta0: f64, psip0: f64, rho0: f64, zeta0: f64, mu: f64) -> Self {
        Self {
            t0,
            theta0,
            psip0,
            rho0,
            zeta0,
            mu,
        }
    }
}
