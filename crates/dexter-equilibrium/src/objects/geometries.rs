//! Representation of an equilibrium's general geometry.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    debug_assert_non_negative_r, equilibrium_type_getter_impl, fluxes_values_array_getter_impl,
    fortran_vec_to_carray2d_impl, lcfs_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl, shape2d_getter_impl,
};
use core::f64::consts::PI;
use dexter_common::array1D_getter_impl;
use ndarray::{Array1, Array2, Order::ColumnMajor};
use rsl_interpolation::{
    Accelerator, Accelerator2d, DynInterpolator, DynInterpolator2d, Interpolation,
    Interpolation1dType, Interpolation2d, Interpolation2dType,
};
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{EqError, EvalError};
use crate::{EquilibriumType, FluxCommute, Geometry};

// ===============================================================================================

/// Analytical Large Aspect Ratio Geometry of a circular device.
///
/// The definitions are not very strict at the moment.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub struct LarGeometry {
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` in [m].
    raxis: f64,
    /// The `r` coordinate's value at the last closed flux surface `rlast` in [m].
    ///
    /// In LAR configuration, `rlast` coincides with the device's minor radius.
    rlast: f64,
    /// The value of the last closed toroidal flux surface `ψ_last` in Normalized units.
    ///
    /// In LAR configuration, the wall coincides with the last closed flux surface.
    psi_last: f64,
}

impl LarGeometry {
    /// Creates a new `LarCurrent`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let geometry = LarGeometry::new(2.0, 1.75, 0.5);
    /// ```
    #[must_use]
    pub fn new(baxis: f64, raxis: f64, rlast: f64) -> Self {
        let psi_last_si = baxis * rlast.powi(2) / 2.0;
        let psi_last = psi_last_si / (baxis * raxis.powi(2));
        Self {
            equilibrium_type: EquilibriumType::Analytical,
            baxis,
            raxis,
            rlast,
            psi_last,
        }
    }
}

impl Geometry for LarGeometry {
    fn psi_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Bad
    }

    fn r_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!((2.0 * psi).sqrt()))
    }

    fn r_of_psip(&self, _: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "r(ψp) (defined through q)".into(),
        ))
    }

    fn psi_of_r(&self, r: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_r!(r);
        Ok(debug_assert_is_finite!(r.powi(2) / 2.0))
    }

    fn psip_of_r(&self, _: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "ψp(r) (defined through q)".into(),
        ))
    }

    fn rlab_of_psi(&self, psi: f64, theta: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(
            self.raxis + self.raxis * (2.0 * psi).sqrt() * theta.cos()
        ))
    }

    fn rlab_of_psip(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "R(ψp, θ) (defined through q)".into(),
        ))
    }

    fn zlab_of_psi(&self, psi: f64, theta: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(
            self.raxis * (2.0 * psi).sqrt() * theta.sin()
        ))
    }

    fn zlab_of_psip(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "Z(ψp, θ) (defined through q)".into(),
        ))
    }

    fn jacobian_of_psi(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "J(ψ, θ) (defined through q, g, I and B)".into(),
        ))
    }

    fn jacobian_of_psip(&self, _: f64, _: f64, _: &mut Accelerator2d) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "J(ψp, θ) (defined through q, g, I and B)".into(),
        ))
    }

    fn rlab_last(&self) -> Array1<f64> {
        let arr = Array1::linspace(0.0, 2.0 * PI, 1000);
        let acc = &mut Accelerator2d::new();
        arr.mapv(|theta| match self.rlab_of_psi(self.psi_last, theta, acc) {
            Ok(rlab_lcfs_value) => rlab_lcfs_value,
            Err(_) => unreachable!("Expression is analytical, cannot fail"),
        })
    }

    fn zlab_last(&self) -> Array1<f64> {
        let arr = Array1::linspace(0.0, 2.0 * PI, 1000);
        let acc = &mut Accelerator2d::new();
        arr.mapv(|theta| match self.zlab_of_psi(self.psi_last, theta, acc) {
            Ok(zlab_lcfs_value) => zlab_lcfs_value,
            Err(_) => unreachable!("Expression is analytical, cannot fail"),
        })
    }
}

impl LarGeometry {
    equilibrium_type_getter_impl!();

    /// Returns the magnetic field strength on the axis `B0` in **\[T\]**.
    #[must_use]
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Returns the horizontal position of the magnetic axis `R0` in **\[m\]**.
    #[must_use]
    pub fn raxis(&self) -> f64 {
        self.raxis
    }

    /// Returns the vertical position of the magnetic axis **in \[m\]**.
    #[must_use]
    pub fn zaxis(&self) -> f64 {
        0.0
    }

    /// Returns the geometrical axis (device major radius) **in \[m\]**.
    #[must_use]
    pub fn rgeo(&self) -> f64 {
        self.raxis
    }

    /// Returns the `r` coordinate's value at the last closed flux surface **in \[m\]**.
    #[must_use]
    pub fn rlast(&self) -> f64 {
        self.rlast
    }

    /// Returns the value of the last closed toroidal flux surface `ψ_last`.
    #[must_use]
    pub fn psi_last(&self) -> f64 {
        self.psi_last
    }
}

// ===============================================================================================

/// Used to create an [`NcGeometry`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
#[derive(Debug)]
pub struct NcGeometryBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// The 1D interpolation type.
    interp1d_type: Interpolation1dType,
    /// The 2D interpolation type.
    interp2d_type: Interpolation2dType,
}

impl NcGeometryBuilder {
    /// Creates a new [`NcGeometryBuilder`] from a netCDF file at `path`, with 1D interpolation
    /// type `interp1d_type` and 2D interpolation type `interp2d_type`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let interp1d_type = Interpolation1dType::Akima;
    /// let interp2d_type = Interpolation2dType::Bicubic;
    /// let builder = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type);
    /// ```
    #[must_use]
    pub fn new(
        path: &Path,
        interp1d_type: Interpolation1dType,
        interp2d_type: Interpolation2dType,
    ) -> Self {
        Self {
            path: path.to_path_buf(),
            interp1d_type,
            interp2d_type,
        }
    }

    /// Creates a new [`NcGeometry`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let interp1d_type = Interpolation1dType::Akima;
    /// let interp2d_type = Interpolation2dType::Bicubic;
    /// let builder = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type);
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EqError`] if it fails to build the [`NcGeometry`].
    pub fn build(self) -> Result<NcGeometry, EqError> {
        NcGeometry::build(self)
    }
}

// ===============================================================================================

/// Describes the general geometry of the equilibrium.
///
/// Stores fluxes, angles and lab variables' data, and provides interpolation methods between them.
///
/// Should be created with an [`NcGeometryBuilder`].
#[non_exhaustive]
#[derive(Clone)]
pub struct NcGeometry {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,

    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// The 1D interpolation type.
    interp1d_type: Interpolation1dType,
    /// The 2D interpolation type.
    interp2d_type: Interpolation2dType,

    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` in [m].
    raxis: f64,
    /// The vertical position of the magnetic axis in [m].
    zaxis: f64,
    /// The geometrical axis (device major radius) in [m].
    rgeo: f64,

    /// The boozer toroidal angle `θ` in [rads].
    theta_values: Vec<f64>,
    /// The toroidal flux coordinate.
    psi: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// The `ψp(ψ)` interpolator.
    psip_of_psi_interp: Option<DynInterpolator>,
    /// The `ψ(ψp)` interpolator.
    psi_of_psip_interp: Option<DynInterpolator>,

    /// The radial coordinate r in [m].
    r_values: Vec<f64>,
    /// The `r(ψ)` interpolator.
    r_of_psi_interp: Option<DynInterpolator>,
    /// The `r(ψp)` interpolator.
    r_of_psip_interp: Option<DynInterpolator>,
    /// The `ψ(r)` interpolator.
    psi_of_r_interp: Option<DynInterpolator>,
    /// The `ψp(r)` interpolator.
    psip_of_r_interp: Option<DynInterpolator>,

    /// The `R` coordinate in [m], flattened in F order.
    rlab_values_fortran_flat: Vec<f64>,
    /// `R(ψ, θ)` interpolator.
    rlab_of_psi_interp: Option<DynInterpolator2d>,
    /// `R(ψp, θ)` interpolator.
    rlab_of_psip_interp: Option<DynInterpolator2d>,

    /// The `Z` coordinate in [m], flattened in F order.
    zlab_values_fortran_flat: Vec<f64>,
    /// `Z(ψ, θ)` interpolator.
    zlab_of_psi_interp: Option<DynInterpolator2d>,
    /// `Z(ψp, θ)` interpolator.
    zlab_of_psip_interp: Option<DynInterpolator2d>,

    /// The VMEC output to Boozer Jacobian in [m/T], flattened in F order.
    jacobian_values_fortran_flat: Vec<f64>,
    /// `J(ψ, θ)` interpolator.
    jacobian_of_psi_interp: Option<DynInterpolator2d>,
    /// `J(ψp, θ)` interpolator.
    jacobian_of_psip_interp: Option<DynInterpolator2d>,
}

/// Creation.
impl NcGeometry {
    /// Constructs an [`NcGeometry`] from an [`NcGeometryBuilder`].
    pub(crate) fn build(builder: NcGeometryBuilder) -> Result<Self, EqError> {
        use crate::extract;
        use crate::extract::netcdf_fields::{NC_BAXIS, NC_RAXIS, NC_RGEO, NC_ZAXIS};
        use crate::extract::netcdf_fields::{NC_JACOBIAN, NC_R, NC_RLAB, NC_THETA, NC_ZLAB};

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let theta_values = extract::array_1d(&file, NC_THETA)?.to_vec();
        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);

        let baxis: f64 = extract::scalar(&file, NC_BAXIS)?;
        let raxis: f64 = extract::scalar(&file, NC_RAXIS)?;
        let zaxis: f64 = extract::scalar(&file, NC_ZAXIS)?;
        let rgeo: f64 = extract::scalar(&file, NC_RGEO)?;
        let r_values = extract::array_1d(&file, NC_R)?.to_vec();
        let rlab_values_fortran_flat = extract::array_2d(&file, NC_RLAB)?
            .flatten_with_order(ColumnMajor)
            .to_vec();
        let zlab_values_fortran_flat = extract::array_2d(&file, NC_ZLAB)?
            .flatten_with_order(ColumnMajor)
            .to_vec();
        let jacobian_values_fortran_flat = extract::array_2d(&file, NC_JACOBIAN)?
            .flatten_with_order(ColumnMajor)
            .to_vec();

        debug_assert!(baxis.is_finite(), "NaN 'baxis' field encountered");
        debug_assert!(raxis.is_finite(), "NaN 'raxis' field encountered");
        debug_assert!(zaxis.is_finite(), "NaN 'zaxis' field encountered");
        debug_assert!(rgeo.is_finite(), "NaN 'rgeo' field encountered");
        debug_assert_all_finite_values(&theta_values);
        debug_assert_all_finite_values(&r_values);
        debug_assert_all_finite_values(&rlab_values_fortran_flat);
        debug_assert_all_finite_values(&zlab_values_fortran_flat);
        debug_assert_all_finite_values(&jacobian_values_fortran_flat);

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let psip_of_psi_interp =
            if (psi.state() == Good) & (psip.state() != FluxCoordinateState::NoValues) {
                Some(DynInterpolator::build(
                    builder.interp1d_type,
                    psi.uvalues(),
                    psip.uvalues(),
                )?)
            } else {
                None
            };
        let psi_of_psip_interp =
            if (psip.state() == Good) & (psi.state() != FluxCoordinateState::NoValues) {
                Some(DynInterpolator::build(
                    builder.interp1d_type,
                    psip.uvalues(),
                    psi.uvalues(),
                )?)
            } else {
                None
            };

        let r_of_psi_interp = match psi.state() {
            Good => DynInterpolator::build(builder.interp1d_type, psi.uvalues(), &r_values).ok(),
            _ => None,
        };
        let r_of_psip_interp = match psip.state() {
            Good => DynInterpolator::build(builder.interp1d_type, psip.uvalues(), &r_values).ok(),
            _ => None,
        };

        // Neither the fluxes or `r` is guaranteed to exist.
        // If `r` exists, then it is guaranteed it's in increasing order.
        let psi_of_r_interp = match psi.state() {
            FluxCoordinateState::NoValues => None,
            _ => DynInterpolator::build(builder.interp1d_type, &r_values, psi.uvalues()).ok(),
        };
        let psip_of_r_interp = match psip.state() {
            FluxCoordinateState::NoValues => None,
            _ => DynInterpolator::build(builder.interp1d_type, &r_values, psip.uvalues()).ok(),
        };

        let rlab_of_psi_interp = match psi.state() {
            Good => DynInterpolator2d::build(
                builder.interp2d_type,
                psi.uvalues(),
                &theta_values,
                &rlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };
        let rlab_of_psip_interp = match psip.state() {
            Good => DynInterpolator2d::build(
                builder.interp2d_type,
                psip.uvalues(),
                &theta_values,
                &rlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };

        let zlab_of_psi_interp = match psi.state() {
            Good => DynInterpolator2d::build(
                builder.interp2d_type,
                psi.uvalues(),
                &theta_values,
                &zlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };
        let zlab_of_psip_interp = match psip.state() {
            Good => DynInterpolator2d::build(
                builder.interp2d_type,
                psip.uvalues(),
                &theta_values,
                &zlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };

        let jacobian_of_psi_interp = match psi.state() {
            Good => DynInterpolator2d::build(
                builder.interp2d_type,
                psi.uvalues(),
                &theta_values,
                &jacobian_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };
        let jacobian_of_psip_interp = match psip.state() {
            Good => DynInterpolator2d::build(
                builder.interp2d_type,
                psip.uvalues(),
                &theta_values,
                &jacobian_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp1d_type: builder.interp1d_type,
            interp2d_type: builder.interp2d_type,
            theta_values,
            psi,
            psip,
            psi_of_psip_interp,
            psip_of_psi_interp,
            baxis,
            raxis,
            zaxis,
            rgeo,
            r_values,
            r_of_psi_interp,
            r_of_psip_interp,
            psi_of_r_interp,
            psip_of_r_interp,
            rlab_values_fortran_flat,
            rlab_of_psi_interp,
            rlab_of_psip_interp,
            zlab_values_fortran_flat,
            zlab_of_psi_interp,
            zlab_of_psip_interp,
            jacobian_values_fortran_flat,
            jacobian_of_psi_interp,
            jacobian_of_psip_interp,
        })
    }
}

impl FluxCommute for NcGeometry {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.psip_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                self.psip.uvalues(),
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψp(ψ)".into())),
        }
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.psi_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                self.psi.uvalues(),
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψ(ψp)".into())),
        }
    }
}

impl Geometry for NcGeometry {
    fn psi_state(&self) -> FluxCoordinateState {
        self.psi.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip.state()
    }

    fn r_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.r_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.r_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("r(ψ)".into())),
        }
    }

    fn r_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.r_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.r_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("r(ψp)".into())),
        }
    }

    fn psi_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_r!(r);
        match self.psi_of_r_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                &self.r_values,
                self.psi.uvalues(),
                r,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψ(r)".into())),
        }
    }

    fn psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_r!(r);
        match self.psip_of_r_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                &self.r_values,
                self.psip.uvalues(),
                r,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψp(r)".into())),
        }
    }

    fn rlab_of_psi(&self, psi: f64, theta: f64, acc: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.rlab_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.rlab_values_fortran_flat,
                psi,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("R(ψ, θ)".into())),
        }
    }

    fn rlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.rlab_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.rlab_values_fortran_flat,
                psip,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("R(ψp, θ)".into())),
        }
    }

    fn zlab_of_psi(&self, psi: f64, theta: f64, acc: &mut Accelerator2d) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.zlab_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.zlab_values_fortran_flat,
                psi,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("Z(ψ, θ)".into())),
        }
    }

    fn zlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.zlab_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.zlab_values_fortran_flat,
                psip,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("Z(ψp, θ)".into())),
        }
    }

    fn jacobian_of_psi(
        &self,
        psi: f64,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.jacobian_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.jacobian_values_fortran_flat,
                psi,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("J(ψ, θ)".into())),
        }
    }

    fn jacobian_of_psip(
        &self,
        psip: f64,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.jacobian_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.jacobian_values_fortran_flat,
                psip,
                theta,
                acc,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("J(ψp, θ)".into())),
        }
    }

    fn rlab_last(&self) -> Array1<f64> {
        // Get the last row of the C-ordered `rlab_array`
        self.rlab_array().row(self.shape().0 - 1).to_owned()
    }

    fn zlab_last(&self) -> Array1<f64> {
        // Get the last row of the C-ordered `zlab_array`
        self.zlab_array().row(self.shape().0 - 1).to_owned()
    }
}

/// Getters.
impl NcGeometry {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();

    /// Returns the 1D interpolation type.
    #[must_use]
    pub fn interp1d_type(&self) -> Interpolation1dType {
        self.interp1d_type
    }

    /// Returns the 2D interpolation type.
    #[must_use]
    pub fn interp2d_type(&self) -> Interpolation2dType {
        self.interp2d_type
    }

    /// Returns the magnetic field strength on the axis `B0` **in \[T\]**.
    #[must_use]
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Returns the horizontal position of the magnetic axis `R0` **in \[m\]**.
    #[must_use]
    pub fn raxis(&self) -> f64 {
        self.raxis
    }

    /// Returns the vertical position of the magnetic axis **in \[m\]**.
    #[must_use]
    pub fn zaxis(&self) -> f64 {
        self.zaxis
    }

    /// Returns the geometrical axis (device major radius) **in \[m\]**.
    #[must_use]
    pub fn rgeo(&self) -> f64 {
        self.rgeo
    }

    /// Returns the `r` coordinate's value at the last closed flux surface **in \[m\]**.
    #[must_use]
    pub fn rlast(&self) -> f64 {
        match self.r_values.last().copied() {
            Some(rlast) => rlast,
            None => unreachable!("NcGeometry cannot be created if `r_values` dont exist"),
        }
    }

    shape2d_getter_impl!();
    lcfs_getter_impl!();
    fluxes_values_array_getter_impl!();
    array1D_getter_impl!(theta_array, theta_values, theta);
    array1D_getter_impl!(r_array, r_values, r);
    fortran_vec_to_carray2d_impl!(rlab_array, rlab_values_fortran_flat, R);
    fortran_vec_to_carray2d_impl!(zlab_array, zlab_values_fortran_flat, Z);
    fortran_vec_to_carray2d_impl!(jacobian_array, jacobian_values_fortran_flat, J);
}

impl std::fmt::Debug for NcGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcGeometry")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("1D interpolation type", &self.interp1d_type())
            .field("2D interpolation type", &self.interp2d_type())
            .field("baxis [T]", &self.baxis)
            .field("raxis [m]", &self.raxis)
            .field("zaxis [m]", &self.zaxis)
            .field("rgeo [m]", &self.rgeo)
            .field("rlast [m]", &self.rlast())
            .field("shape (ψ/ψp, θ)", &self.shape())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_geometry(path_str: &str) -> NcGeometry {
        let path = PathBuf::from(&path_str);
        let builder = NcGeometryBuilder::new(
            &path,
            Interpolation1dType::Steffen,
            Interpolation2dType::Bicubic,
        );
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
        let geometry = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(geometry.psi_state(), FluxCoordinateState::Good);
        assert_eq!(geometry.psip_state(), FluxCoordinateState::Bad);

        assert_eq!(geometry.psi.state(), FluxCoordinateState::Good);
        assert_eq!(geometry.psip.state(), FluxCoordinateState::Bad);
        assert!(geometry.psip_of_psi_interp.is_some());
        assert!(geometry.r_of_psi_interp.is_some());
        assert!(geometry.rlab_of_psi_interp.is_some());
        assert!(geometry.zlab_of_psi_interp.is_some());
        assert!(geometry.jacobian_of_psi_interp.is_some());

        assert!(geometry.psi_of_psip_interp.is_none());
        assert!(geometry.r_of_psip_interp.is_none());
        assert!(geometry.rlab_of_psip_interp.is_none());
        assert!(geometry.zlab_of_psip_interp.is_none());
        assert!(geometry.jacobian_of_psip_interp.is_none());

        assert!(geometry.psi_of_r_interp.is_some());
        assert!(geometry.psip_of_r_interp.is_some());

        assert!(geometry.psi_array().is_some());
        assert!(geometry.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psi_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        assert!(g.psip_of_psi(p, acc1).unwrap().is_finite());
        assert!(g.r_of_psi(p, acc1).unwrap().is_finite());
        assert!(g.rlab_of_psi(p, t, acc2).unwrap().is_finite());
        assert!(g.zlab_of_psi(p, t, acc2).unwrap().is_finite());
        assert!(g.jacobian_of_psi(p, t, acc2).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        use EvalError::UndefinedEvaluation as err;
        matches!(g.psi_of_psip(p, acc1), Err(err(..)));
        matches!(g.r_of_psip(p, acc1), Err(err(..)));
        matches!(g.rlab_of_psip(p, t, acc2), Err(err(..)));
        matches!(g.zlab_of_psip(p, t, acc2), Err(err(..)));
        matches!(g.jacobian_of_psip(p, t, acc2), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let geometry = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(geometry.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(geometry.psip_state(), FluxCoordinateState::Good);

        assert_eq!(geometry.psi.state(), FluxCoordinateState::Bad);
        assert_eq!(geometry.psip.state(), FluxCoordinateState::Good);
        assert!(geometry.psip_of_psi_interp.is_none());
        assert!(geometry.r_of_psi_interp.is_none());
        assert!(geometry.rlab_of_psi_interp.is_none());
        assert!(geometry.zlab_of_psi_interp.is_none());
        assert!(geometry.jacobian_of_psi_interp.is_none());

        assert!(geometry.psi_of_psip_interp.is_some());
        assert!(geometry.r_of_psip_interp.is_some());
        assert!(geometry.rlab_of_psip_interp.is_some());
        assert!(geometry.zlab_of_psip_interp.is_some());
        assert!(geometry.jacobian_of_psip_interp.is_some());

        assert!(geometry.psi_of_r_interp.is_some());
        assert!(geometry.psip_of_r_interp.is_some());

        assert!(geometry.psi_array().is_some());
        assert!(geometry.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psip_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        assert!(g.psi_of_psip(p, acc1).unwrap().is_finite());
        assert!(g.r_of_psip(p, acc1).unwrap().is_finite());
        assert!(g.rlab_of_psip(p, t, acc2).unwrap().is_finite());
        assert!(g.zlab_of_psip(p, t, acc2).unwrap().is_finite());
        assert!(g.jacobian_of_psip(p, t, acc2).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let p = 0.01;
        let t = 3.14;
        use EvalError::UndefinedEvaluation as err;
        matches!(g.psip_of_psi(p, acc1), Err(err(..)));
        matches!(g.r_of_psi(p, acc1), Err(err(..)));
        matches!(g.rlab_of_psi(p, t, acc2), Err(err(..)));
        matches!(g.zlab_of_psi(p, t, acc2), Err(err(..)));
        matches!(g.jacobian_of_psi(p, t, acc2), Err(err(..)));
    }
}
