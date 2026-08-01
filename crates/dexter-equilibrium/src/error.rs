//! Custom Error types for netCDF and equilibrium object creation/evalution errors.

/// Top level Error type.
#[derive(thiserror::Error, Debug)]
pub enum EqError {
    /// From [`std::io::Error`].
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// From [`NcError`].
    #[error("NetCDF error: {0}")]
    NcError(#[from] NcError),

    /// From [`EvalError`].
    #[error("Eval error: {0}")]
    EvalError(#[from] EvalError),

    /// Error from [`rsl_interpolation`].
    #[error("Interpolation error: {0}")]
    InterpolationError(#[from] rsl_interpolation::InterpolatorError),

    /// Analytical threshold index cannot be greater that the number of data points.
    #[error("Analytical threshold index cannot be greater that the number of data points")]
    InvalidFluteModeAnalyticalThresholdIndex,
}

/// Evaluation related errors.
#[derive(thiserror::Error, Debug)]
pub enum EvalError {
    /// 1D Interpolation domain error from [`rsl_interpolation`].
    #[error("1D Interpolation domain error: {0}")]
    Domain1dError(#[from] rsl_interpolation::Domain1dError),

    /// 2D Interpolation domain error from [`rsl_interpolation`].
    #[error("2D Interpolation domain error: {0}")]
    Domain2dError(#[from] rsl_interpolation::Domain2dError),

    /// Analytical evaluation method received an out-of-bounds input.
    ///
    /// The bounds check is enforced by the definition of the equilibrium object's
    /// [`LastClosedFluxSurface`](crate::LastClosedFluxSurface), rather than the formula itself.
    #[error("Analytical calculation domain error")]
    AnalyticalDomainError,

    /// Called undefined evaluation method.
    #[error("Call to undefined evaluation method '{0}'")]
    UndefinedEvaluation(Box<str>),
}

/// netCDF handling errors.
#[derive(thiserror::Error, Debug)]
pub enum NcError {
    /// File not found.
    #[error("File '{0}' not found")]
    FileNotFound(std::path::PathBuf),

    /// Error opening [`netcdf::File`].
    #[error("Error opening '{path}' ({err})")]
    FileOpenError {
        /// Path to file.
        path: std::path::PathBuf,
        /// The wrapped [`netcdf::Error`].
        err: netcdf::Error,
    },

    /// Attribute not found.
    #[error("Attribute '{0}' not found")]
    AttributeNotFound(Box<str>),

    /// Error extracting the value of an Attribute.
    #[error("Error extracting the value of the Attribute {0}")]
    AttributeValueError(Box<str>),

    /// Error extracting netCDF Convention Version.
    #[error("Error extracting netCDF convention version: {0}")]
    VersionError(#[from] semver::Error),

    /// Error extracting values from a [`netcdf::Variable`].
    #[error("Error extracting values from variable {name}: {err}")]
    GetValues {
        /// The name of the variable.
        name: String,
        /// The wrapped [`netcdf::Error`].
        err: netcdf::Error,
    },

    /// Shape mismatch while extracting value into an array.
    #[error("Shape mismatch while extracting value into an array: {0}")]
    ShapeError(#[from] ndarray::ShapeError),

    /// Empty [`netcdf::Variable`].
    #[error("'{0}' variable is empty")]
    EmptyVariable(String),

    /// [`netcdf::Variable`] not found in [`netcdf::File`].
    #[error("'{0}' variable not found in NetCDF file")]
    VariableNotFound(Box<str>),

    /// Mode with passed mode number does not exist.
    #[error("Mode '{which}={mode}' not found in the NetCDF file.")]
    FluteModeNotFound {
        /// The name of the mode ('m' or 'n').
        which: String,
        /// The mode number.
        mode: i64,
    },

    /// θ-padding error.
    #[error("Padding error: {0}")]
    PaddingError(Box<str>),
}
