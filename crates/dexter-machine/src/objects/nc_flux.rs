//! Definition of the `NcFlux` object.

use ndarray::{Array1, ArrayRef1, Axis};

use crate::extract;
use crate::extract::netcdf_fields::{NC_PSI_NORM, NC_PSIP_NORM};

/// Defines whether or not a Flux coordinate can be used for evaluations.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[expect(clippy::exhaustive_enums, reason = "No other possible states")]
pub enum FluxCoordinateState {
    /// Good coordinate: values exist and are increasing. Can be used as an x coordinate for
    /// evaluations.
    Good,
    /// Bad coordinate: values exist but are not monotonic. Can only be used as y-data.
    Bad,
    /// Does not exist or is empty (fallback to when an extraction error occurs).
    NoValues,
}

impl FluxCoordinateState {
    /// Returns `Good` if the values are strictly monotonic, and `Bad`
    /// otherwise.
    #[expect(clippy::float_cmp, reason = "signum() always returns exactly 1.0_f64")]
    fn from_array(arr: &ArrayRef1<f64>) -> Self {
        // Use abs() to handle negative q-factors; we only care for monotonicity
        if arr
            .abs()
            .diff(1, Axis(0))
            .iter()
            .all(|diff| diff.signum() == 1.0_f64)
        {
            Self::Good
        } else {
            Self::Bad
        }
    }
}

/// Representation of a Flux Coordinate created from a netCDF file. Contains the Flux's values
/// and state.
#[derive(Clone)]
pub(crate) struct NcFlux {
    /// The extracted flux values, if the extraction was successful.
    values: Option<Box<[f64]>>,
    /// The state of the coordinate.
    state: FluxCoordinateState,
}

impl NcFlux {
    /// Creates the toroidal flux `ψ` representation.
    pub(crate) fn toroidal(file: &netcdf::File) -> Self {
        Self::build(file, NC_PSI_NORM)
    }

    /// Creates the poroidal flux `ψp` representation.
    pub(crate) fn poloidal(file: &netcdf::File) -> Self {
        Self::build(file, NC_PSIP_NORM)
    }

    /// Creates the Flux Coordinate from raw values.
    pub(crate) fn from_raw_values(values: &[f64]) -> Self {
        Self {
            values: Some(values.to_vec().into_boxed_slice()),
            state: FluxCoordinateState::from_array(&Array1::from(values.to_vec())),
        }
    }

    /// Extracts the values and sets the state.
    ///
    /// It is not guaranteed that both fluxes (good or bad) exist in a dataset, so we dont want to
    /// return an error in case the corresponding field is not found in the netCDF file.
    fn build(file: &netcdf::File, netcdf_field: &str) -> Self {
        match extract::array_1d(file, netcdf_field) {
            Ok(array) => Self {
                values: Some(array.to_vec().into_boxed_slice()),
                state: FluxCoordinateState::from_array(&array),
            },
            Err(_) => Self {
                values: None,
                state: FluxCoordinateState::NoValues,
            },
        }
    }

    /// Returns the state of the flux coordinate.
    pub(crate) fn state(&self) -> FluxCoordinateState {
        self.state
    }

    /// Returns the flux values, if the exist.
    pub(crate) fn values(&self) -> Option<&[f64]> {
        self.values.as_deref()
    }

    /// Returns a slice to the wrapped values, panicking if they do not exist.
    ///
    /// Since the values are referenced by the interpolators all the time, we do the unwrapping
    /// here to avoid messing up the interpolation calls, which results in very ugly code.
    ///
    /// This method should only be used in places where we know the values exists, for example
    /// under an `state != FluxCoordinateState::NoValues` guard.
    pub(crate) fn uvalues(&self) -> &[f64] {
        match self.values.as_ref() {
            Some(values) => values,
            None => unreachable!("Should only be called when values are guaranteed to exist"),
        }
    }

    /// Returns the flux value at the last closed flux surface.
    pub(crate) fn last_value(&self) -> Option<f64> {
        self.values
            .as_ref()
            .and_then(|values| values.last().copied())
    }
}

impl std::fmt::Debug for NcFlux {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let len = match self.values.as_ref() {
            Some(values) => format!("{}", values.len()),
            None => String::from("No values"),
        };
        let last = match self.last_value() {
            Some(value) => format!("{value}"),
            None => String::from("No value"),
        };
        f.debug_struct("Flux")
            .field("state", &self.state)
            .field("len", &len)
            .field("last value", &last)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::extract::{TEST_NETCDF_PATH, open};

    use super::*;

    fn create_empty_flux() -> NcFlux {
        NcFlux {
            values: None,
            state: FluxCoordinateState::NoValues,
        }
    }

    #[test]
    #[should_panic]
    fn get_uvalues_empty() {
        let flux = create_empty_flux();
        let _ = format!("{flux:?}");
        let _ = flux.uvalues();
    }

    #[test]
    fn get_values_empty() {
        let flux = create_empty_flux();
        let _ = format!("{flux:?}");
        assert!(flux.values().is_none());
    }

    #[test]
    fn non_existing_flux_variable() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let file = open(&path).unwrap();
        let noflux = NcFlux::build(&file, "NOT_A_FIELD");
        assert_eq!(noflux.state, FluxCoordinateState::NoValues);
        assert!(noflux.values.is_none());
    }
}
