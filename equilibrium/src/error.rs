use std::path::PathBuf;

/// Custom error types from equilibrium objects
#[derive(thiserror::Error, Debug)]
pub enum EqError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// From [`ndarray::ShapeError`].
    #[error("ShapeError error: {0}")]
    ShapeError(#[from] ndarray::ShapeError),

    /// Error from [`rsl_interpolation`].
    #[error("Interpolation error: {0}")]
    InterpolationError(#[from] rsl_interpolation::InterpolationError),

    /// Interpolation domain error from [`rsl_interpolation`].
    #[error("Interpolation domain error: {0}")]
    DomainError(#[from] rsl_interpolation::DomainError),

    /// From [`NcError`].
    #[error("NetCDF error: {0}")]
    NcError(#[from] NcError),
}

/// Custom error types from extraction methods
#[derive(thiserror::Error, Debug)]
pub enum NcError {
    /// File not found.
    #[error("File '{0}' not found")]
    FileNotFound(PathBuf),

    /// Error opening NetCDF [`netcdf::File`].
    #[error("Error opening '{path}' ({err})")]
    FileOpenError { path: PathBuf, err: netcdf::Error },

    /// Error extracting values from a [`netcdf::Variable`].
    #[error("Error extracting values from variable {name}: {err}")]
    GetValues { name: String, err: netcdf::Error },

    /// Shape mismatch while extracting value into an array.
    #[error("Shape mismatch while extracting value into an array: {0}")]
    NdArray(#[from] ndarray::ShapeError),

    /// Empty [`netcdf::Variable`].
    #[error("'{0}' variable is empty")]
    EmptyVariable(String),

    /// [`netcdf::Variable`] not found in [`netcdf::File`].
    #[error("'{0}' variable not found in NetCDF file")]
    VariableNotFound(Box<str>),

    /// Harmonic with passed m, n mode numbers does not exist.
    #[error("Mode '{which}={mode}' not found in the NetCDF file.")]
    HarmonicModeNotFound { which: String, mode: i64 },
}
