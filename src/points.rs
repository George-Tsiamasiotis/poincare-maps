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
    /// The canonical momentum `Pζ`.
    pub pzeta: f64,
}

#[derive(Debug)]
#[allow(dead_code)]
pub(crate) struct Point {
    /// The time of evaluation.
    pub(crate) t: f64,
    /// The `θ` angle.
    pub(crate) theta: f64,
    /// The poloidal magnetic flux `ψ_p`.
    pub(crate) psip: f64,
    /// The parallel gyro radius `ρ`.
    pub(crate) rho: f64,
    /// The `ζ` angle.
    pub(crate) zeta: f64,
}
