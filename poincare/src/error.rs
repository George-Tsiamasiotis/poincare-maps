#[derive(thiserror::Error, Debug)]
pub enum PoincareError {
    #[error("Particle Error: {0}")]
    ParticleError(#[from] particle::ParticleError),

    #[error("Initial conditions arrays must be of the same size")]
    InitMismatch,
}
