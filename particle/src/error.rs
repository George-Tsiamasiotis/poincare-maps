use std::time::Duration;

#[derive(thiserror::Error, Debug)]
pub enum ParticleError {
    /// Error from [`equilibrium`].
    #[error("Equilibrium error: {0}")]
    EqError(#[from] equilibrium::EqError),

    /// Error from [`rsl_interpolation`].
    #[error("Interpolation error: {0}")]
    InterpolationError(#[from] rsl_interpolation::InterpolationError),

    /// Interpolation domain error.
    #[error("Interpolation domain error: {0}")]
    DomainError(#[from] rsl_interpolation::DomainError),

    /// Particle timed out.
    #[error("Particle timed out after {0:?}")]
    TimedOut(Duration),

    /// Intersection accuracy check failed.
    #[error("Intersection accuracy check failed.")]
    IntersectionError,
}

impl From<ParticleError> for Box<str> {
    fn from(value: ParticleError) -> Self {
        value.to_string().into_boxed_str()
    }
}
