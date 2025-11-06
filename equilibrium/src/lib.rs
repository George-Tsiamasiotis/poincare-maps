mod bfield;
mod currents;
mod error;
mod harmonic;
mod perturbation;
mod qfactor;

pub use bfield::Bfield;
pub use currents::Currents;
pub use error::EqError;
pub use harmonic::{Harmonic, HarmonicCache};
pub use perturbation::Perturbation;
pub use qfactor::Qfactor;

pub type Result<T> = std::result::Result<T, EqError>;

/// Magnetic flux, in Normalized Units.
#[doc(alias = "f64")]
pub type Flux = f64;

/// Angle in radians.
#[doc(alias = "f64")]
pub type Radians = f64;
