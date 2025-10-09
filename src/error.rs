#[derive(thiserror::Error, Debug)]
pub enum MapError {
    /// Error Initializing System.
    #[error("Error Initializing System")]
    SystemInitError,

    /// Error from [`rsl_interpolation`].
    #[error("Interpolation Error: {0}")]
    InterpolationError(#[from] rsl_interpolation::InterpolationError),

    /// Interpolation domain error.
    #[error("Interpolation domain error: {0}")]
    DomainError(#[from] rsl_interpolation::DomainError),

    /// Error from [`tokamak_netcdf`].
    #[error("netCDF error: {0}")]
    NcError(#[from] tokamak_netcdf::NcError),

    /// Supplied angle must be either "zeta" or "theta".
    #[error("Supplied angle must be either 'zeta' or 'theta'")]
    InvalidAngle,
}
