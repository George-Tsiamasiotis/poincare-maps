mod error;
mod initials;

pub use error::PoincareError;
pub use initials::PoincareInit;

pub type Result<T> = std::result::Result<T, PoincareError>;
