pub(crate) type PointTuple = (f64, f64, f64, f64, f64);

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Point {
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

impl Point {
    pub(crate) fn to_tuple(&self) -> PointTuple {
        (self.t, self.theta, self.psip, self.rho, self.zeta)
    }
}
