#[derive(thiserror::Error, Debug)]
pub enum EqError {
    /// Error from [`tokamak_netcdf`].
    #[error("netCDF error: {0}")]
    NcError(#[from] tokamak_netcdf::NcError),

    /// /// Error Initializing System.
    /// #[error("Error Initializing System")]
    /// SystemInitError,

    /// Error from [`rsl_interpolation`].
    #[error("Interpolation Error: {0}")]
    InterpolationError(#[from] rsl_interpolation::InterpolationError),

    /// Interpolation domain error.
    #[error("Interpolation domain error: {0}")]
    DomainError(#[from] rsl_interpolation::DomainError),
}
