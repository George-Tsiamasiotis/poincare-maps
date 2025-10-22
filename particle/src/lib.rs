mod config;
mod error;
mod initial;
mod particle;
mod solver;
mod state;

pub use config::get_config;
pub use error::ParticleError;
pub use initial::InitialConditions;
pub use particle::Particle;
pub use state::State;

pub type Result<T> = std::result::Result<T, ParticleError>;
