mod error;
mod init;

pub use error::PoincareError;
pub use init::PoincareInit;

pub type Result<T> = std::result::Result<T, PoincareError>;
