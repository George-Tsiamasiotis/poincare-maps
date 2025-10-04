mod equilibrium;
mod error;
mod points;
mod rk45;
mod state;
mod system;

pub use equilibrium::Bfield;
pub use equilibrium::Current;
pub use equilibrium::Qfactor;

pub use error::MapError;
pub use points::InitialConditions;
pub use system::System;

pub(crate) use points::Point;
pub(crate) use rk45::Rk45State;
pub(crate) use state::State;

pub type Result<T> = std::result::Result<T, MapError>;
