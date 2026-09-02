//! Representation of an equilibrium's magnetic field.

use crate::{
    EquilibriumObject, debug_assert_is_2pi_modulo, debug_assert_is_finite,
    debug_assert_non_negative_psi, debug_assert_non_negative_psip, fluxes_values_array_getter_impl,
    interp_type_getter_impl, lcfs_getter_impl, netcdf_path_getter_impl, netcdf_version_getter_impl,
};
use ndarray::{Array1, Array2, Axis, Order::ColumnMajor};
use ndarray::{concatenate, s};
use rsl_interpolation::Accelerator2d;
use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::Interpolation2dType;
use crate::constants::DEFAULT_THETA_PADDING_WIDTH;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{Bfield, ObjectType};
use crate::{EqError, EvalError, NcError};
use dexter_common::{DynInterpolator2d, make_interp2d};

// ===============================================================================================

/// Analytical Large Aspect Ratio magnetic field with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).
///
/// The LAR magnetic field can only be expressed with respect to the toroidal flux ψ.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
#[derive(Clone)]
pub struct LarBfield;

impl LarBfield {
    /// Creates a new `LarBfield`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let bfield = LarBfield::new();
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl EquilibriumObject for LarBfield {
    fn object_type(&self) -> ObjectType {
        ObjectType::Analytical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Bad
    }
}

impl Bfield for LarBfield {
    fn b_of_psi(&self, psi: f64, theta: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(
            1.0 - (2.0 * psi).sqrt() * theta.cos()
        ))
    }

    fn b_of_psip(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation("B(ψp, θ)".into()))
    }

    fn db_dpsi(&self, psi: f64, theta: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(-theta.cos() / (2.0 * psi).sqrt()))
    }

    fn db_dpsip(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation("dB(ψp, θ)/dψp".into()))
    }

    fn db_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        _: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!((2.0 * psi).sqrt() * theta.sin()))
    }

    fn db_of_psip_dtheta(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation("dB(ψp, θ)/dθ".into()))
    }
}

impl std::fmt::Debug for LarBfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Large Aspect Ratio Bfield with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).")
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcBfield`].
#[non_exhaustive]
#[derive(Debug)]
pub struct NcBfieldBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// The interpolation type.
    interp_type: Interpolation2dType,
    /// The number of columns to pad the `B` array.
    padding: usize,
}

impl NcBfieldBuilder {
    /// Creates a new [`NcBfieldBuilder`] from a netCDF file at `path`, with 2D interpolation
    /// type `interp_type`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic);
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: Interpolation2dType) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type,
            padding: DEFAULT_THETA_PADDING_WIDTH,
        }
    }

    /// Sets the left-right `θ` padding width.
    ///
    /// At the grid edges, the interpolator's higher derivatives are not well defined. By
    /// left-right padding the `B` array with extra `ψ=const` columns, we force the interpolator
    /// to take `θ`'s periodicity into account and therefore calculate the correct derivative
    /// values.
    ///
    /// Note that in contrast to the one-dimensional cubic spline, in a bicubic interpolation 3
    /// columns are not enough to ensure periodicity, since the spline coefficients depend on the
    /// values of the whole array.
    ///
    /// According to [`this`] stack overflow thread, the effect of the `i`-th column at the `j`-th
    /// column of the spline scales as `r^|i-j|`, where `r = sqrt(3)-2 ≈ -0.26`. Therefore, with a
    /// padding of 10, the effect at the `θ=0` boundary would be of the order of 1e-6.
    ///
    /// # Default
    ///
    /// The default padding value is 15 columns, where the relative error is about `1e-9`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let typ = Interpolation2dType::Bicubic;
    /// let builder = NcBfieldBuilder::new(&path, typ).with_padding(5).build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// [`this`]: https://stackoverflow.com/a/25106574/32596387
    #[must_use]
    pub fn with_padding(mut self, padding: usize) -> Self {
        self.padding = padding;
        self
    }

    /// Creates a new [`NcBfield`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let typ = Interpolation2dType::Bicubic;
    /// let bfield = NcBfieldBuilder::new(&path, typ).build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EqError`] if it fails to build the [`NcBfield`].
    pub fn build(self) -> Result<NcBfield, EqError> {
        NcBfield::build(self)
    }
}

// ===============================================================================================

/// Numerical magnetic field profile reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcBfieldBuilder`].
#[non_exhaustive]
pub struct NcBfield {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,
    /// The interpolation type.
    interp_type: Interpolation2dType,

    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The number of columns to pad the `B` array.
    padding: usize,

    /// The boozer toroidal angle `θ` in [rads], as extracted from the netCDF file.
    theta_values: Box<[f64]>,
    /// The boozer toroidal angle `θ` in [rads], with the added padding.
    theta_values_padded: Box<[f64]>,
    /// The toroidal flux coordinate.
    psi: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// The `B` array as extracted from the netCDF file.
    b_array: Array2<f64>,
    /// The `B` values, flattened in F order, with the added padding.
    b_values_fortran_flat_padded: Box<[f64]>,
    /// `B(ψ, θ)` interpolator.
    b_of_psi_interp: Option<DynInterpolator2d>,
    /// `B(ψp, θ)` interpolator.
    b_of_psip_interp: Option<DynInterpolator2d>,
}

/// Creation.
impl NcBfield {
    /// Constructs an [`NcBfield`] from an [`NcBfieldBuilder`].
    pub(crate) fn build(builder: NcBfieldBuilder) -> Result<Self, EqError> {
        use crate::extract;
        use crate::extract::netcdf_fields::{NC_B_NORM, NC_BAXIS, NC_THETA};

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);

        let theta_array = extract::array_1d(&file, NC_THETA)?;
        let b_array = extract::array_2d(&file, NC_B_NORM)?;
        let baxis: f64 = extract::scalar(&file, NC_BAXIS)?;

        let theta_values_padded = Self::pad_theta_array(&theta_array, builder.padding)?.to_vec();
        let b_array_padded = Self::pad_b_array(&b_array, builder.padding)?;
        let b_values_fortran_flat_padded = b_array_padded.flatten_with_order(ColumnMajor).to_vec();

        debug_assert!(baxis.is_finite(), "NaN 'baxis' field encountered");
        debug_assert_all_finite_values(&theta_values_padded);
        debug_assert_all_finite_values(&b_values_fortran_flat_padded);

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let b_of_psi_interp = match psi.state() {
            Good => Some(make_interp2d(
                builder.interp_type,
                psi.uvalues(),
                &theta_values_padded,
                &b_values_fortran_flat_padded,
            )?),
            _ => None,
        };
        let b_of_psip_interp = match psip.state() {
            Good => Some(make_interp2d(
                builder.interp_type,
                psip.uvalues(),
                &theta_values_padded,
                &b_values_fortran_flat_padded,
            )?),
            _ => None,
        };

        Ok(Self {
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            theta_values: theta_array.to_vec().into_boxed_slice(),
            theta_values_padded: theta_values_padded.into_boxed_slice(),
            psi,
            psip,
            baxis,
            padding: builder.padding,
            b_array,
            b_values_fortran_flat_padded: b_values_fortran_flat_padded.into_boxed_slice(),
            b_of_psi_interp,
            b_of_psip_interp,
        })
    }

    /// Returns the left-right padded `θ` array.
    fn pad_theta_array(theta_array: &Array1<f64>, padding: usize) -> Result<Array1<f64>, NcError> {
        if padding > theta_array.len() {
            return Err(NcError::PaddingError(
                "'padding' cannot be bigger that the number of 'θ' values.".into(),
            ));
        }

        let left_values = theta_array.slice(s![1..=padding]).to_owned();
        let right_values = theta_array
            .slice(s![(-1 - padding as isize)..=-2])
            .to_owned();
        let left_pad = right_values - TAU;
        let right_pad = left_values + TAU;

        concatenate(
            Axis(0),
            &[left_pad.view(), theta_array.view(), right_pad.view()],
        )
        .or(Err(NcError::PaddingError("Concatenation error.".into())))
    }

    /// Returns the left-right padded `B` array.
    fn pad_b_array(b_array: &Array2<f64>, padding: usize) -> Result<Array2<f64>, NcError> {
        if padding > b_array.ncols() {
            return Err(NcError::PaddingError(
                "'padding' cannot be bigger that the number of 'θ' values.".into(),
            ));
        }
        let right_padding = b_array.slice(s![.., 1..=padding]);
        let left_padding = b_array.slice(s![.., (-1 - padding as isize)..=-2]);

        concatenate(
            Axis(1),
            &[left_padding.view(), b_array.view(), right_padding.view()],
        )
        .or(Err(NcError::PaddingError("Concatenation error.".into())))
    }
}

impl EquilibriumObject for NcBfield {
    fn object_type(&self) -> ObjectType {
        ObjectType::Numerical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        self.psi.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip.state()
    }
}

impl Bfield for NcBfield {
    fn b_of_psi(&self, psi: f64, theta: f64, acc: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        debug_assert_is_2pi_modulo!(theta);
        match self.b_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values_padded,
                &self.b_values_fortran_flat_padded,
                psi,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("B(ψ, θ)".into())),
        }
    }

    fn b_of_psip(&self, psip: f64, theta: f64, acc: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        debug_assert_is_2pi_modulo!(theta);
        match self.b_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values_padded,
                &self.b_values_fortran_flat_padded,
                psip,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("B(ψp, θ)".into())),
        }
    }

    fn db_dpsi(&self, psi: f64, theta: f64, acc: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        debug_assert_is_2pi_modulo!(theta);
        match self.b_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv_x(
                self.psi.uvalues(),
                &self.theta_values_padded,
                &self.b_values_fortran_flat_padded,
                psi,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dB(ψ, θ)/dψ".into())),
        }
    }

    fn db_dpsip(&self, psip: f64, theta: f64, acc: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        debug_assert_is_2pi_modulo!(theta);
        match self.b_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv_x(
                self.psip.uvalues(),
                &self.theta_values_padded,
                &self.b_values_fortran_flat_padded,
                psip,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dB(ψp, θ)/dψp".into())),
        }
    }

    fn db_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        debug_assert_is_2pi_modulo!(theta);
        match self.b_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv_y(
                self.psi.uvalues(),
                &self.theta_values_padded,
                &self.b_values_fortran_flat_padded,
                psi,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dB(ψ, θ)/dθ".into())),
        }
    }

    fn db_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        debug_assert_is_2pi_modulo!(theta);
        match self.b_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv_y(
                self.psip.uvalues(),
                &self.theta_values_padded,
                &self.b_values_fortran_flat_padded,
                psip,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dB(ψp, θ)/dθ".into())),
        }
    }
}

/// Getters.
impl NcBfield {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    interp_type_getter_impl!(2);

    /// Returns the magnetic field strength on the axis `B0` **in \[T\]**.
    #[must_use]
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Returns the number of `θ` padding columns (per side).
    #[must_use]
    pub fn padding(&self) -> usize {
        self.padding
    }

    /// Returns the `θ` array, as extracted from the netCDF file.
    #[must_use]
    pub fn theta_array(&self) -> Array1<f64> {
        Array1::from(self.theta_values.clone())
    }

    /// Returns the padded `θ` array.
    #[must_use]
    pub fn theta_array_padded(&self) -> Array1<f64> {
        Array1::from(self.theta_values_padded.clone())
    }

    /// Returns the `B` array, as extracted from the netCDF file.
    #[must_use]
    pub fn b_array(&self) -> Array2<f64> {
        self.b_array.clone()
    }

    /// Returns the padded `B` array.
    #[must_use]
    pub fn b_array_padded(&self) -> Array2<f64> {
        // Array is in Fortran order, so we must reverse the shape
        let shape = (
            self.b_array.ncols() + 2 * self.padding,
            self.b_array.nrows(),
        );
        #[expect(clippy::missing_panics_doc, reason = "infallible")]
        Array2::from_shape_vec(shape, self.b_values_fortran_flat_padded.to_vec())
            .expect("Shape is correct by definition")
            .reversed_axes()
    }

    /// Returns the (ψ/ψp, θ) shape of the **padded** arrays that were used to create the interpolator.
    #[must_use]
    pub fn shape_padded(&self) -> (usize, usize) {
        let mut actual_shape = self.shape();
        actual_shape.1 += 2 * self.padding;
        actual_shape
    }

    /// Returns the (ψ/ψp, θ) shape of the initial 1D arrays (before the padding).
    #[must_use]
    pub fn shape(&self) -> (usize, usize) {
        let psi_len = match self.psi_state() {
            FluxCoordinateState::NoValues => 0,
            _ => self.psi.uvalues().len(),
        };
        let psip_len = match self.psip_state() {
            FluxCoordinateState::NoValues => 0,
            _ => self.psip.uvalues().len(),
        };
        let xlen = psi_len.max(psip_len);
        (xlen, self.theta_values.len())
    }

    lcfs_getter_impl!();
    fluxes_values_array_getter_impl!();
}

impl std::fmt::Debug for NcBfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcBfield")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("interpolation type", &self.interp_type())
            .field("baxis [T]", &self.baxis)
            .field("shape (ψ/ψp, θ)", &self.shape())
            .field("padding", &self.padding)
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_bfield(path_str: &str) -> NcBfield {
        let path = PathBuf::from(&path_str);
        let builder = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic);
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
        let bfield = create_nc_bfield(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(bfield.psi_state(), FluxCoordinateState::Good);
        assert_eq!(bfield.psip_state(), FluxCoordinateState::Bad);

        assert_eq!(bfield.psi_state(), FluxCoordinateState::Good);
        assert_eq!(bfield.psip_state(), FluxCoordinateState::Bad);
        assert!(bfield.b_of_psi_interp.is_some());
        assert!(bfield.b_of_psip_interp.is_none());

        assert!(bfield.psi_array().is_some());
        assert!(bfield.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let b = create_nc_bfield(TOROIDAL_TEST_NETCDF_PATH);
        let acc = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        assert!(b.b_of_psi(p, t, acc).unwrap().is_finite());
        assert!(b.db_dpsi(p, t, acc).unwrap().is_finite());
        assert!(b.db_of_psi_dtheta(p, t, acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let b = create_nc_bfield(TOROIDAL_TEST_NETCDF_PATH);
        let acc = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        use EvalError::UndefinedEvaluation as err;
        matches!(b.b_of_psip(p, t, acc), Err(err(..)));
        matches!(b.db_dpsip(p, t, acc), Err(err(..)));
        matches!(b.db_of_psip_dtheta(p, t, acc), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let bfield = create_nc_bfield(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(bfield.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(bfield.psip_state(), FluxCoordinateState::Good);

        assert_eq!(bfield.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(bfield.psip_state(), FluxCoordinateState::Good);
        assert!(bfield.b_of_psi_interp.is_none());
        assert!(bfield.b_of_psip_interp.is_some());

        assert!(bfield.psi_array().is_some());
        assert!(bfield.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let b = create_nc_bfield(POLOIDAL_TEST_NETCDF_PATH);
        let acc = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        assert!(b.b_of_psip(p, t, acc).unwrap().is_finite());
        assert!(b.db_dpsip(p, t, acc).unwrap().is_finite());
        assert!(b.db_of_psip_dtheta(p, t, acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let b = create_nc_bfield(POLOIDAL_TEST_NETCDF_PATH);
        let acc = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        use EvalError::UndefinedEvaluation as err;
        matches!(b.b_of_psi(p, t, acc), Err(err(..)));
        matches!(b.db_dpsi(p, t, acc), Err(err(..)));
        matches!(b.db_of_psi_dtheta(p, t, acc), Err(err(..)));
    }
}

#[cfg(test)]
mod lar_values {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    #[rustfmt::skip]
    fn lar_bfield_values_gcmotion_check() {
        let b = LarBfield::new();
        let acc = &mut Accelerator2d::new();
        let epsilon = 1e-20;

        let (psi, theta) = (0.2, 0.4);
        assert_relative_eq!(b.b_of_psi(psi, theta, acc).unwrap(), 0.4174698790024389, epsilon=epsilon);
        assert_relative_eq!(b.db_dpsi(psi, theta, acc).unwrap(), -1.4563253024939025, epsilon=epsilon);
        assert_relative_eq!(b.db_of_psi_dtheta(psi, theta, acc).unwrap(), 0.24628978486848968, epsilon=epsilon);

        let (psi, theta) = (15.0, 1000.0);
        assert_relative_eq!(b.b_of_psi(psi, theta, acc).unwrap(), -2.0802770595333673, epsilon=epsilon);
        assert_relative_eq!(b.db_dpsi(psi, theta, acc).unwrap(), -0.10267590198444558, epsilon=epsilon);
        assert_relative_eq!(b.db_of_psi_dtheta(psi, theta, acc).unwrap(), 4.529005766888851, epsilon=epsilon);
    }
}

#[cfg(test)]
mod padding {
    use ndarray::array;

    use super::*;

    #[test]
    #[rustfmt::skip]
    fn theta_padding() {
        let theta_values = Array1::from(vec![0.0, 2.0, 4.0, TAU]);
        let pad0 = NcBfield::pad_theta_array(&theta_values, 0).unwrap().to_vec();
        let pad1 = NcBfield::pad_theta_array(&theta_values, 1).unwrap().to_vec();
        let pad2 = NcBfield::pad_theta_array(&theta_values, 2).unwrap().to_vec();
        let pad3 = NcBfield::pad_theta_array(&theta_values, 3).unwrap().to_vec();

        assert_eq!(pad0, theta_values.to_vec());
        assert_eq!(pad1, vec![
                4.0-TAU,
                0.0, 2.0, 4.0, TAU,
                2.0+TAU
            ]
        );
        assert_eq!(pad2, vec![
                2.0-TAU, 4.0-TAU,
                0.0, 2.0, 4.0, TAU,
                2.0+TAU, 4.0+TAU
            ]
        );
        assert_eq!(pad3, vec![
                -TAU, 2.0-TAU, 4.0-TAU,
                0.0, 2.0, 4.0, TAU,
                2.0+TAU, 4.0+TAU, TAU+TAU
            ]
        );
    }

    #[test]
    fn b_padding() {
        let b_array = array![
            [1.0, 4.0, 7.0, 1.0],
            [2.0, 5.0, 8.0, 2.0],
            [3.0, 6.0, 9.0, 3.0]
        ];
        let pad0 = NcBfield::pad_b_array(&b_array, 0).unwrap();
        let pad1 = NcBfield::pad_b_array(&b_array, 1).unwrap();
        let pad2 = NcBfield::pad_b_array(&b_array, 2).unwrap();

        assert_eq!(pad0, b_array);
        assert_eq!(
            pad1,
            array![
                [7.0, 1.0, 4.0, 7.0, 1.0, 4.0],
                [8.0, 2.0, 5.0, 8.0, 2.0, 5.0],
                [9.0, 3.0, 6.0, 9.0, 3.0, 6.0]
            ]
        );
        assert_eq!(
            pad2,
            array![
                [4.0, 7.0, 1.0, 4.0, 7.0, 1.0, 4.0, 7.0],
                [5.0, 8.0, 2.0, 5.0, 8.0, 2.0, 5.0, 8.0],
                [6.0, 9.0, 3.0, 6.0, 9.0, 3.0, 6.0, 9.0]
            ]
        );
    }
}
