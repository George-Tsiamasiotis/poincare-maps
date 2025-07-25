#[derive(thiserror::Error)]
///Custom Error types
pub enum SplineError {
    /// One of the supplied datasets is empty.
    #[error("Supplied datasets cannot be empty.")]
    InvalidDataset,

    /// `x` points dataset is not sorted.
    #[error("Supplied x dataset must be sorted.")]
    UnsortedDataset,

    /// 'x' and `y` datasets have differnet length.
    #[error("Supplied datasets must be 1D and of equal length.")]
    DatasetMismatch,

    /// Supplied array size is less than the interpolation type's minimum size.
    #[error("Supplied array size is less than the interpolation type's minimum size.")]
    NotEnoughPoints,

    /// Error calling `gsl_interp_accel_alloc`.
    #[error("Error calling gsl_interp_accel_alloc.")]
    GSLAccelAlloc,

    /// Error calling `gsl_interp_alloc`. `rgsl`'s return an Option, so no error.
    #[error("Error calling gsl_interp_alloc.")]
    GSLInterpAlloc,

    /// Error calling 'gsl_interp_init'.
    #[error("Error calling gsl_interp_init: {err:?}.")]
    GSLSplineInit { err: rgsl::Value },

    /// Supplied x is out of bounds. GSL crashes hard when this happens so we catch it earlier.
    #[error("Supplied x out of bounds (GSL error: {err:?}).")]
    GSLInputDomainError { err: rgsl::Value },
}

impl std::fmt::Debug for SplineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{self}")
    }
}
