pub(crate) use std::{path::PathBuf, time::Duration};

#[derive(thiserror::Error, Debug)]
pub enum MapError {
    /// Cannot file netCDF file at `path`
    #[error("Cannot find netCDF file at `{path}`: {source}")]
    PathError {
        #[source]
        source: std::io::Error,
        path: PathBuf,
    },

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

    /// Error running particle.
    #[error("Error running particle: {0}")]
    OrbitError(Box<str>),

    /// Particle integration time out.
    #[error("Particle timed out after {0:?} and {1} steps.")]
    OrbitTimeout(Duration, usize),
}
