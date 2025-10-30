use crate::Point;

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
    /// Creates a new set of initial conditions.
    ///
    /// # Example
    ///
    /// ```
    /// # use particle::*;
    /// #
    /// # fn main() {
    /// # let initial = InitialConditions{
    ///     t0: 0.0,
    ///     theta0: 3.14,
    ///     psip0: 0.015,
    ///     rho0: 0.001,
    ///     zeta0: 0.0,
    ///     mu: 0.0,
    /// };
    /// # }
    /// ```
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

    /// Creates a [`Point`] from `self`.
    ///
    /// The derived quantites Pζ, Pθ, and ψ are uninitialized.
    ///
    /// # Example
    ///
    /// ```
    /// # use particle::*;
    /// #
    /// # fn main() {
    /// # let initial = InitialConditions{
    ///     t0: 0.0,
    ///     theta0: 3.14,
    ///     psip0: 0.015,
    ///     rho0: 0.001,
    ///     zeta0: 0.0,
    ///     mu: 0.0,
    /// };
    /// let point: Point = initial.to_point();
    /// # }
    /// ```
    pub fn to_point(&self) -> Point {
        Point {
            time: self.t0,
            theta: self.theta0,
            psip: self.psip0,
            rho: self.rho0,
            zeta: self.zeta0,
            mu: self.mu,
            ..Default::default()
        }
    }
}
