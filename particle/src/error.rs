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

    /// TODO:
    #[error("TODO")]
    InvalidAngle,
}
