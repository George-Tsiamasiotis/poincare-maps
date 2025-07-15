//! This crate provides a cleaner and more suitable [`NetCDF`] interface for use with reconstructed equilibria
//! from [`Tokamak`] devices.
//!
//!
//! ## Example
//!
//! ```no_run
//! # use std::path::PathBuf;
//! # use tokamak_netcdf::NcError;
//! #
//! # fn main() -> Result<(), NcError> {
//! #
//!     // Path must be relative to the directory where "cargo run" is called
//!     let path = PathBuf::from(r"./reconstructed/data.nc");
//!     let nc_data = tokamak_netcdf::NcData::open(path)?;
//!
//!     println!("{:#?}", nc_data);
//!
//! # Ok(())
//! # }
//! ```
//!
//! ## Note
//!
//! This crate requires the [`netCDF-C`] library, which is available in most linux package managers. In
//! case it is not, it can be statically linked with the 'static' feature, which is provided by the
//! [`netcdf crate`].
//!
//! [`NetCDF`]: https://www.unidata.ucar.edu/software/netcdf
//! [`netCDF-C`]: https://github.com/Unidata/netcdf-c
//! [`netcdf crate`]: https://github.com/georust/netcdf
//! [`Tokamak`]: https://en.wikipedia.org/wiki/Tokamak

mod error;
mod extract;
mod open;

mod bfield;
mod coords;
mod currents;
mod scalars;

pub use error::NcError;
pub use open::NcData;

pub type Result<T> = std::result::Result<T, NcError>;
