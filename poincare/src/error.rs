#[derive(thiserror::Error, Debug)]
pub enum PoincareError {
    #[error("{0}")]
    ParticleError(#[from] particle::ParticleError),

    #[error("Initial conditions arrays must be of the same size")]
    InitMismatch,
}
