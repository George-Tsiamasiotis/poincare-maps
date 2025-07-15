use std::path::PathBuf;

#[derive(thiserror::Error)]
/// The error type for opening and extracting data from NetCDF files.
pub enum NcError {
    /// Catches errors from the wrapped netcdf library.
    /// Lower level error: netcdf::Error::Netcdf.
    #[error("Wrapped library Error: {reason} ({source}).")]
    LibraryError {
        #[source]
        source: netcdf::Error,
        reason: Box<str>,
    },

    /// Given path is not a file.
    #[error("'{0}': File not found.")]
    FileNotFound(PathBuf),

    /// Variable does not exist.
    #[error("'{0}' variable not found.")]
    VariableNotFound(Box<str>),

    /// Variable exists, but has no data.
    #[error("'{0}' variable is empty.")]
    EmptyVariable(Box<str>),

    /// Expected scalar value, found array.
    #[error("'{0}' variable is not scalar.")]
    NotScalar(Box<str>),

    /// Attempted to extract non-1D variable as 1D.
    #[error("'{0}' variable is not 1-dimensional.")]
    Not1D(Box<str>),

    /// Attempted to extract non-2D variable as 2D.
    #[error("'{0}' variable is not 2-dimensional.")]
    Not2D(Box<str>),

    /// Errors from netcdf's `get_values()` functions. Those are hard to track but should be
    /// basically unreachable.
    #[error("Error extracting values from '{name}' variable: {source}.")]
    GetValuesError {
        #[source]
        source: netcdf::Error,
        name: Box<str>,
    },

    /// Supplied wrong `Extents` dimensionality to `get_values()`.
    /// Lower level error: netcdf::Error::DimensionalityMismatch.
    #[error("Error extracting '{name}': {source}.")]
    DimensionError {
        #[source]
        source: netcdf::Error,
        name: Box<str>,
    },
}

impl std::fmt::Debug for NcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{self}")
    }
}
