/// Custom error types
#[derive(thiserror::Error, Debug)]
pub enum EqError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Error from [`ndarray`].
    #[error("ShapeError error: {0}")]
    ShapeError(#[from] ndarray::ShapeError),

    /// Error from [`tokamak_netcdf`].
    #[error("NetCDF error: {0}")]
    NcError(#[from] tokamak_netcdf::NcError),

    /// Error from [`rsl_interpolation`].
    #[error("Interpolation error: {0}")]
    InterpolationError(#[from] rsl_interpolation::InterpolationError),

    /// Interpolation domain error from [`rsl_interpolation`]..
    #[error("Interpolation domain error: {0}")]
    DomainError(#[from] rsl_interpolation::DomainError),
}
