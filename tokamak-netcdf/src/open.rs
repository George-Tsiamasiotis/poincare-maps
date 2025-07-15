//! Handles NetCDF file opening and `NcData` creation.

use std::path::PathBuf;

use crate::bfield::Bfield;
use crate::coords::Coords;
use crate::currents::Currents;
use crate::scalars::Scalars;
use crate::{NcError, Result};

/// NetCDF equilibrium data.
pub struct NcData {
    /// Path to NetCDF file.
    pub path: PathBuf,
    /// Equilibrium's scalar values.
    pub scalars: Scalars,
    /// Equilibrium's coordinate variables.
    pub coords: Coords,
    /// Representation of the plasma toroidal **I** and poloidal **g** currents.
    pub currents: Currents,
    /// Representation of the equilibrium's magnetic field **B**.
    pub bfield: Bfield,
}

impl NcData {
    /// Creates an NcData from a NetCDF file.
    ///
    /// ## Example
    ///
    /// ```no_run
    /// # use std::path::PathBuf;
    /// # use tokamak_netcdf::NcError;
    /// #
    /// # fn main() -> Result<(), NcError> {
    /// #
    /// // Path must be relative to the directory where "cargo run" is called
    /// let path = PathBuf::from(r"./reconstructed/data.nc");
    /// let nc_data = tokamak_netcdf::NcData::open(path)?;
    ///
    /// println!("{:#?}", nc_data);
    /// # Ok(())
    /// # }
    /// ```
    pub fn open(path: PathBuf) -> Result<Self> {
        use NcError::*;

        if !path.exists() {
            return Err(FileNotFound(path));
        }

        // If this fails, its due to an underlying library error.
        let nc_file = match netcdf::open(&path) {
            Ok(nc_file) => nc_file,
            Err(liberror) => {
                return Err(LibraryError {
                    source: liberror, // Error::Netcdf
                    reason: "Error opening NetCDF file".into(),
                });
            }
        };

        let scalars = Scalars::build(&nc_file)?;
        let coords = Coords::build(&nc_file)?;
        let currents = Currents::build(&nc_file)?;
        let bfield = Bfield::build(&nc_file)?;

        let rec = NcData {
            path,
            scalars,
            coords,
            currents,
            bfield,
        };

        Ok(rec)
    }
}

impl std::fmt::Debug for NcData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcData")
            .field("scalars", &self.scalars)
            .field("coords", &self.coords)
            .field("currents", &self.currents)
            .field("bfield", &self.bfield)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use crate::NcData;
    use crate::NcError;

    #[test]
    fn file_not_found() {
        let x: Result<NcData, NcError> = NcData::open("not_an_existing_path".into());
        assert!(matches!(x.unwrap_err(), NcError::FileNotFound { .. }));
    }
}
