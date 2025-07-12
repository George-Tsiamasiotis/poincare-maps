mod acc;
mod error;
mod evals;
mod interp_types;
mod spline;

// Types refering to rgsl's Interp and Spline types, to avoid confusion.
pub(crate) type RgslSpline = rgsl::Spline;
pub(crate) type RgslInterpType = rgsl::InterpType;
pub(crate) type RgslInterpAccel = rgsl::InterpAccel;

pub use acc::Accelerator;
pub use error::SplineError;
pub use interp_types::InterpolationType;
pub use spline::Spline;

pub type Result<T, E = SplineError> = std::result::Result<T, E>;
