mod error;
mod evolution;
mod initial;
mod mapping;
mod particle;
mod point;
mod rkf45;
mod state;

pub use error::ParticleError;
pub use evolution::Evolution;
pub use initial::InitialConditions;
pub use mapping::*;
pub use particle::{IntegrationStatus, Particle};
pub use point::Point;
pub use state::State;

pub(crate) use rkf45::Solver;

pub type Result<T> = std::result::Result<T, ParticleError>;
