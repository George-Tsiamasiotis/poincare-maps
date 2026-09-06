//! Representation of a machine's plasma current.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_flux, fluxes_values_array_getter_impl,
    interp_type_getter_impl, lcfs_getter_impl, netcdf_path_getter_impl, netcdf_version_getter_impl,
};
use core::hint::cold_path;
use ndarray::Array1;
use rsl_interpolation::Accelerator;
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::Interpolation1dType;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{Current, MachineObject, MachineType};
use crate::{EvalError, MachineError};
use crate::{MagneticFlux, MagneticFlux::*};
use dexter_common::{DynInterpolator, array1D_getter_impl, make_interp};

// ===============================================================================================

/// Analytical Large Aspect Ratio Current with g=1 and I=0.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
#[derive(Clone)]
pub struct LarCurrent;

impl LarCurrent {
    /// Creates a new `LarCurrent`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let current = LarCurrent::new();
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl MachineObject for LarCurrent {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }
}

impl Current for LarCurrent {
    fn eval_g(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(1.0)
    }

    fn eval_i(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(0.0)
    }

    fn eval_g_deriv(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(0.0)
    }

    fn eval_i_deriv(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(0.0)
    }
}

impl std::fmt::Debug for LarCurrent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Large Aspect Ratio Current with g=1 and I=0")
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcCurrent`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct NcCurrentBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// The interpolation type.
    interp_type: Interpolation1dType,
}

impl NcCurrentBuilder {
    /// Creates a new `NcCurrentBuilder` from a netCDF file at `path`, with `interp_type`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcCurrentBuilder::new(&path, Interpolation1dType::Cubic);
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: Interpolation1dType) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type,
        }
    }

    /// Creates a new [`NcCurrent`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let current = NcCurrentBuilder::new(&path, Interpolation1dType::Akima).build()?;
    /// Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`MachineError`] if it fails to build the [`NcCurrent`].
    pub fn build(self) -> Result<NcCurrent, MachineError> {
        NcCurrent::build(self)
    }
}

// ===============================================================================================

/// Numerical plasma current reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcCurrentBuilder`].
#[non_exhaustive]
pub struct NcCurrent {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,
    /// The interpolation type.
    interp_type: Interpolation1dType,

    /// The toroidal flux coordinate.
    psi: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// The `g` values.
    g_values: Box<[f64]>,
    /// The `g(ψ)` interpolatior.
    g_of_psi_interp: Option<DynInterpolator>,
    /// The `g(ψp)` interpolatior.
    g_of_psip_interp: Option<DynInterpolator>,

    /// `I` values.
    i_values: Box<[f64]>,
    /// The `I(ψ)` interpolatior.
    i_of_psi_interp: Option<DynInterpolator>,
    /// The `I(ψp)` interpolatior.
    i_of_psip_interp: Option<DynInterpolator>,
}

/// Creation.
impl NcCurrent {
    /// Constructs an `NcCurrent` from [`NcCurrentBuilder`].
    fn build(builder: NcCurrentBuilder) -> Result<Self, MachineError> {
        use crate::extract;
        use crate::extract::netcdf_fields::{NC_G_NORM, NC_I_NORM};

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);
        let g_values = extract::array_1d(&file, NC_G_NORM)?.to_vec();
        let i_values = extract::array_1d(&file, NC_I_NORM)?.to_vec();

        debug_assert_all_finite_values(&g_values);
        debug_assert_all_finite_values(&i_values);

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let g_of_psi_interp = match psi.state() {
            Good => Some(make_interp(builder.interp_type, psi.uvalues(), &g_values)?),
            _ => None,
        };
        let i_of_psi_interp = match psi.state() {
            Good => Some(make_interp(builder.interp_type, psi.uvalues(), &i_values)?),
            _ => None,
        };

        let g_of_psip_interp = match psip.state() {
            Good => Some(make_interp(builder.interp_type, psip.uvalues(), &g_values)?),
            _ => None,
        };
        let i_of_psip_interp = match psip.state() {
            Good => Some(make_interp(builder.interp_type, psip.uvalues(), &i_values)?),
            _ => None,
        };

        Ok(Self {
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            psi,
            psip,
            g_values: g_values.into_boxed_slice(),
            i_values: i_values.into_boxed_slice(),
            g_of_psi_interp,
            g_of_psip_interp,
            i_of_psi_interp,
            i_of_psip_interp,
        })
    }
}

impl MachineObject for NcCurrent {
    fn machine_type(&self) -> MachineType {
        MachineType::Numerical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        self.psi.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip.state()
    }
}

impl Current for NcCurrent {
    fn eval_g(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let (x, xa, interp) = match flux {
            Toroidal(v) => (v, self.psi.uvalues(), self.g_of_psi_interp.as_ref()),
            Poloidal(v) => (v, self.psip.uvalues(), self.g_of_psip_interp.as_ref()),
        };
        let ya = &self.g_values;
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(interp.eval(xa, ya, x, acc)?))
        } else {
            cold_path();
            let msg = format!("g({})", flux.kind());
            Err(EvalError::UndefinedEvaluation(msg))
        }
    }

    fn eval_i(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let (x, xa, interp) = match flux {
            Toroidal(v) => (v, self.psi.uvalues(), self.i_of_psi_interp.as_ref()),
            Poloidal(v) => (v, self.psip.uvalues(), self.i_of_psip_interp.as_ref()),
        };
        let ya = &self.i_values;
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(interp.eval(xa, ya, x, acc)?))
        } else {
            cold_path();
            let msg = format!("I({})", flux.kind());
            Err(EvalError::UndefinedEvaluation(msg))
        }
    }

    fn eval_g_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let (x, xa, interp) = match flux {
            Toroidal(v) => (v, self.psi.uvalues(), self.g_of_psi_interp.as_ref()),
            Poloidal(v) => (v, self.psip.uvalues(), self.g_of_psip_interp.as_ref()),
        };
        let ya = &self.g_values;
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(interp.eval_deriv(xa, ya, x, acc)?))
        } else {
            cold_path();
            let msg = format!("dg({})/d{}", flux.kind(), flux.kind());
            Err(EvalError::UndefinedEvaluation(msg))
        }
    }

    fn eval_i_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let (x, xa, interp) = match flux {
            Toroidal(v) => (v, self.psi.uvalues(), self.i_of_psi_interp.as_ref()),
            Poloidal(v) => (v, self.psip.uvalues(), self.i_of_psip_interp.as_ref()),
        };
        let ya = &self.i_values;
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(interp.eval_deriv(xa, ya, x, acc)?))
        } else {
            cold_path();
            let msg = format!("dI({})/d{}", flux.kind(), flux.kind());
            Err(EvalError::UndefinedEvaluation(msg))
        }
    }
}

/// Getters.
impl NcCurrent {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    interp_type_getter_impl!(1);
    lcfs_getter_impl!();
    fluxes_values_array_getter_impl!();
    array1D_getter_impl!(g_array, g_values, g);
    array1D_getter_impl!(i_array, i_values, I);
}

impl std::fmt::Debug for NcCurrent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcCurrent")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("interpolation type", &self.interp_type())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_current(path_str: &str) -> NcCurrent {
        let path = PathBuf::from(path_str);
        let builder = NcCurrentBuilder::new(&path, Interpolation1dType::Cubic);
        builder.build().unwrap()
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(current.psi_state(), FluxCoordinateState::Good);
        assert_eq!(current.psip_state(), FluxCoordinateState::Bad);

        assert_eq!(current.psi.state(), FluxCoordinateState::Good);
        assert_eq!(current.psip.state(), FluxCoordinateState::Bad);
        assert!(current.g_of_psi_interp.is_some());
        assert!(current.i_of_psi_interp.is_some());
        assert!(current.g_of_psip_interp.is_none());
        assert!(current.i_of_psip_interp.is_none());

        assert!(current.psi_array().is_some());
        assert!(current.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Toroidal(0.01);
        assert!(current.eval_g(flux, &mut acc).unwrap().is_finite());
        assert!(current.eval_i(flux, &mut acc).unwrap().is_finite());
        assert!(current.eval_i_deriv(flux, &mut acc).unwrap().is_finite());
        assert!(current.eval_i_deriv(flux, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Poloidal(0.01);
        use EvalError::UndefinedEvaluation as err;
        matches!(current.eval_g(flux, &mut acc), Err(err(..)));
        matches!(current.eval_i(flux, &mut acc), Err(err(..)));
        matches!(current.eval_i_deriv(flux, &mut acc), Err(err(..)));
        matches!(current.eval_i_deriv(flux, &mut acc), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(current.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(current.psip_state(), FluxCoordinateState::Good);

        assert_eq!(current.psi.state(), FluxCoordinateState::Bad);
        assert_eq!(current.psip.state(), FluxCoordinateState::Good);
        assert!(current.g_of_psi_interp.is_none());
        assert!(current.i_of_psi_interp.is_none());
        assert!(current.g_of_psip_interp.is_some());
        assert!(current.i_of_psip_interp.is_some());

        assert!(current.psi_array().is_some());
        assert!(current.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Poloidal(0.01);
        assert!(current.eval_g(flux, &mut acc).unwrap().is_finite());
        assert!(current.eval_i(flux, &mut acc).unwrap().is_finite());
        assert!(current.eval_i_deriv(flux, &mut acc).unwrap().is_finite());
        assert!(current.eval_i_deriv(flux, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Toroidal(0.01);
        use EvalError::UndefinedEvaluation as err;
        matches!(current.eval_g(flux, &mut acc), Err(err(..)));
        matches!(current.eval_i(flux, &mut acc), Err(err(..)));
        matches!(current.eval_i_deriv(flux, &mut acc), Err(err(..)));
        matches!(current.eval_i_deriv(flux, &mut acc), Err(err(..)));
    }
}
