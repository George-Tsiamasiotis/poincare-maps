#![doc = include_str!("../README.md")]

mod bfield;
mod currents;
mod error;
mod harmonic;
mod perturbation;
mod qfactor;

pub mod extract;

pub use bfield::Bfield;
pub use currents::Currents;
pub use error::{EqError, NcError};
pub use harmonic::{Harmonic, HarmonicCache};
pub use perturbation::Perturbation;
pub use qfactor::Qfactor;

pub use config::STUB_NETCDF_PATH;
pub use config::netcdf_fields;

pub type Result<T> = std::result::Result<T, EqError>;

/// Magnetic flux, in Normalized Units.
#[doc(alias = "f64")]
pub type Flux = f64;

/// Angle in radians.
#[doc(alias = "f64")]
pub type Radians = f64;

/// Distance, in Normalized Units (normalized to the major radius R).
#[doc(alias = "f64")]
pub type Length = f64;
