#[derive(thiserror::Error, Debug)]
pub enum PoincareError {
    #[error("Initial conditions arrays must be of the same size")]
    InitMismatch,
}
