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
