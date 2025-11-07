/// Custom error types
#[derive(thiserror::Error, Debug)]
pub enum ParticleError {
    /// Error from [`equilibrium`].
    #[error("{0}")]
    EqError(#[from] equilibrium::EqError),

    /// Particle timed out.
    #[error("Particle timed out after {0:?}")]
    TimedOut(std::time::Duration),

    /// Intersection accuracy check failed.
    #[error("Intersection accuracy check failed.")]
    IntersectionError,
}
