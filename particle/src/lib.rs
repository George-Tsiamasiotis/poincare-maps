mod config;
mod error;
mod evolution;
mod particle;
mod point;
mod rkf45;
mod state;

pub use config::get_config;
pub use error::ParticleError;
pub use evolution::Evolution;
pub use particle::{IntegrationStatus, Particle};
pub use point::Point;
pub use state::State;

pub(crate) use rkf45::Solver;

pub type Result<T> = std::result::Result<T, ParticleError>;
