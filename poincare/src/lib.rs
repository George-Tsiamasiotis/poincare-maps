#![doc = include_str!("../README.md")]

mod error;
mod initials;
mod poincare;
mod results;

pub use error::PoincareError;
pub use initials::PoincareInit;
pub use poincare::Poincare;
pub use results::PoincareResults;

pub type Result<T> = std::result::Result<T, PoincareError>;

pub use equilibrium::{Flux, Radians};
pub use particle::{Length, MagneticMoment, Time};
