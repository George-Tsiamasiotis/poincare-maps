mod bfield;
mod current;
mod error;
mod harmonic;
mod perturbation;
mod qfactor;

pub use bfield::Bfield;
pub use current::Current;
pub use error::EqError;
pub use harmonic::Harmonic;
pub use perturbation::Perturbation;
pub use qfactor::Qfactor;

pub type Result<T> = std::result::Result<T, EqError>;
