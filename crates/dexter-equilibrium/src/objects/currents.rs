//! Representation of an equilibrium's plasma current.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    equilibrium_type_getter_impl, fluxes_values_array_getter_impl, interp_type_getter_impl,
    lcfs_getter_impl, netcdf_path_getter_impl, netcdf_version_getter_impl,
};
use dexter_common::array1D_getter_impl;
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolator, Interpolation, Interpolation1dType};
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{Current, EquilibriumType};
use crate::{EqError, EvalError};

// ===============================================================================================

/// Analytical Large Aspect Ratio Current with g=1 and I=0.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
#[derive(Clone)]
pub struct LarCurrent {
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
}

impl LarCurrent {
    /// Creates a new `LarCurrent`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let current = LarCurrent::new();
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self {
            equilibrium_type: EquilibriumType::Analytical,
        }
    }

    equilibrium_type_getter_impl!();
}

impl Current for LarCurrent {
    fn psi_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn g_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(1.0)
    }

    fn g_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(1.0)
    }

    fn i_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(0.0)
    }

    fn i_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(0.0)
    }

    fn dg_dpsi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(0.0)
    }

    fn dg_dpsip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(0.0)
    }

    fn di_dpsi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(0.0)
    }

    fn di_dpsip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
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
    /// Creates a new [`NcCurrentBuilder`] from a netCDF file at `path`, with `interp_type`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
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
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let current = NcCurrentBuilder::new(&path, Interpolation1dType::Akima).build()?;
    /// Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EqError`] if it fails to build the [`NcCurrent`].
    pub fn build(self) -> Result<NcCurrent, EqError> {
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
#[derive(Clone)]
pub struct NcCurrent {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,

    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
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
    /// Constructs an [`NcCurrent`] from [`NcCurrentBuilder`].
    fn build(builder: NcCurrentBuilder) -> Result<Self, EqError> {
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
            Good => Some(DynInterpolator::build(
                builder.interp_type,
                psi.uvalues(),
                &g_values,
            )?),
            _ => None,
        };
        let i_of_psi_interp = match psi.state() {
            Good => Some(DynInterpolator::build(
                builder.interp_type,
                psi.uvalues(),
                &i_values,
            )?),
            _ => None,
        };

        let g_of_psip_interp = match psip.state() {
            Good => Some(DynInterpolator::build(
                builder.interp_type,
                psip.uvalues(),
                &g_values,
            )?),
            _ => None,
        };
        let i_of_psip_interp = match psip.state() {
            Good => Some(DynInterpolator::build(
                builder.interp_type,
                psip.uvalues(),
                &i_values,
            )?),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
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

impl Current for NcCurrent {
    fn psi_state(&self) -> FluxCoordinateState {
        self.psi.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip.state()
    }

    fn g_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.g_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.g_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("g(ψ)".into())),
        }
    }

    fn g_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.g_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.g_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("g(ψp)".into())),
        }
    }

    fn i_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.i_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.i_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("I(ψ)".into())),
        }
    }

    fn i_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.i_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.i_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("I(ψp)".into())),
        }
    }

    fn dg_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.g_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psi.uvalues(),
                &self.g_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dg(ψ)/dψ".into())),
        }
    }

    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.g_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psip.uvalues(),
                &self.g_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dg(ψp)/dψp".into())),
        }
    }

    fn di_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.i_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psi.uvalues(),
                &self.i_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dI(ψ)/dψ".into())),
        }
    }

    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.i_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psip.uvalues(),
                &self.i_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dI(ψp)/dψp".into())),
        }
    }
}

/// Getters.
impl NcCurrent {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
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
            .field("equilibrium type", &self.equilibrium_type())
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
        assert!(current.g_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(current.i_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(current.dg_dpsi(0.01, &mut acc).unwrap().is_finite());
        assert!(current.di_dpsi(0.01, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EvalError::UndefinedEvaluation as err;
        matches!(current.g_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(current.i_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(current.dg_dpsip(0.01, &mut acc), Err(err(..)));
        matches!(current.di_dpsip(0.01, &mut acc), Err(err(..)));
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
        assert!(current.g_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(current.i_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(current.dg_dpsip(0.01, &mut acc).unwrap().is_finite());
        assert!(current.di_dpsip(0.01, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EvalError::UndefinedEvaluation as err;
        matches!(current.g_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(current.i_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(current.dg_dpsi(0.01, &mut acc), Err(err(..)));
        matches!(current.di_dpsi(0.01, &mut acc), Err(err(..)));
    }
}
