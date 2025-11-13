/// Custom error types
#[derive(thiserror::Error, Debug)]
pub enum ParticleError {
    /// From [`equilibrium::EqError`].
    #[error("{0}")]
    EqError(#[from] equilibrium::EqError),

    /// From [`equilibrium::NcError`].
    #[error("{0}")]
    NcError(#[from] equilibrium::NcError),

    /// Particle timed out.
    #[error("Particle timed out after {0:?}")]
    TimedOut(std::time::Duration),

    /// Intersection accuracy check failed.
    #[error("Intersection accuracy check failed.")]
    IntersectionError,

    /// NaN encountered inside solver.
    #[error("NaN encountered inside solver.")]
    SolverNan,
}
