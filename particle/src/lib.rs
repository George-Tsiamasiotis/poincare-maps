#![doc = include_str!("../README.md")]

mod error;
mod evolution;
mod mapping;
mod particle;
mod rkf45;
mod state;

pub use error::ParticleError;
pub use evolution::Evolution;
pub use mapping::*;
pub use particle::{InitialConditions, IntegrationStatus, Particle};
pub use state::State;

pub(crate) use rkf45::Solver;

pub type Result<T> = std::result::Result<T, ParticleError>;

pub use equilibrium::Flux;
pub use equilibrium::Radians;

/// Distance, in Normalized Units (normalized to the major radius R).
#[doc(alias = "f64")]
pub type Length = f64;

/// Time, in Normalized Units (inversed gyrofrequency on magnetic axis).
#[doc(alias = "f64")]
pub type Time = f64;

/// Magnetic Moment, in Normalized Units.
#[doc(alias = "f64")]
pub type MagneticMoment = f64;

/// Canonical Momentum, in Normalized Units.
#[doc(alias = "f64")]
pub type CanonicalMomentum = f64;

/// Energy, in Normalized Units.
#[doc(alias = "f64")]
pub type Energy = f64;
