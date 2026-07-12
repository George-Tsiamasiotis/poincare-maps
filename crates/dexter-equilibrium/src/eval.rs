//! Definitions of evaluation methods of equilibrium objects.
//!
//! For analytical equilibria, this is achieved by evaluation of analytical formulas, while for
//! numerical equilibria by interpolation over the reconstructed data arrays.

use ndarray::Array1;
use rsl_interpolation::{Accelerator, Cache};

use crate::{EvalError, FluxCoordinateState};

/// Reference to a dynamically dispatched [`Mode`] object.
pub type DynMode = Box<dyn Mode>;

/// Reference to a dynamically dispatched [`ModeCache`] object.
pub type DynModeCache = Box<dyn ModeCache>;

/// Equilibrium geometry related quantities computation.
pub trait Geometry {
    /// Returns the [`FluxCoordinateState`] of the toroidal `ψ` flux coordinate.
    fn psi_state(&self) -> FluxCoordinateState;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψp` flux coordinate.
    fn psip_state(&self) -> FluxCoordinateState;

    /// Calculates the radial coordinate `r(ψ)` in **\[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let r_of_psi = geometry.r_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn r_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the radial coordinate `r(ψp)` in **\[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let r_of_psip = geometry.r_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn r_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `ψ(r)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psi_of_r = geometry.psi_of_r(0.02, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn psi_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `ψp(r)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psip_of_r = geometry.psip_of_r(0.02, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `R(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let rlab_of_psi = geometry.rlab_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn rlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `R(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psip_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let rlab_of_psip = geometry.rlab_of_psip(0.01, 3.14, &mut psip_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn rlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `Z(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let zlab_of_psi = geometry.zlab_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn zlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `Z(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psip_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let zlab_of_psip = geometry.zlab_of_psip(0.01, 3.14, &mut psip_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn zlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates the Jacobian `J(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let jacobian_of_psi = geometry.jacobian_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn jacobian_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates the Jacobian `J(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psip_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let jacobian_of_psip = geometry.zlab_of_psip(0.01, 3.14, &mut psip_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn jacobian_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Returns the last `Rlab` values that correspond to the device's last closed flux surface.
    fn rlab_last(&self) -> Array1<f64>;

    /// Returns the last `Zlab` values that correspond to the device's last closed flux surface.
    fn zlab_last(&self) -> Array1<f64>;
}

/// Conversion between the two flux coordinates `ψ` and `ψp`.
pub trait FluxCommute {
    /// Calculates `ψp(ψ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psip_of_psi = qfactor.psip_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the conversion fails for any reason.
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `ψ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psi_of_psip = geometry.psi_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the conversion fails for any reason.
    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;
}

/// q-factor related quantities computation.
pub trait Qfactor {
    /// Returns the [`FluxCoordinateState`] of the toroidal `ψ` flux coordinate.
    fn psi_state(&self) -> FluxCoordinateState;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψp` flux coordinate.
    fn psip_state(&self) -> FluxCoordinateState;

    /// Returns the value of the last closed toroidal flux `ψ_last`.
    fn psi_last(&self) -> f64;

    /// Returns the value of the last closed poloidal flux `ψp_last`.
    fn psip_last(&self) -> f64;

    /// Returns the qfactor's value at the last closed flux surface.
    fn qlast(&self) -> f64;

    /// Returns the qfactor's value at the magnetic axis.
    fn qaxis(&self) -> f64;

    /// Calculates `q(ψ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let q_of_psi = qfactor.q_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `q(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let q_of_psip = qfactor.q_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the derivative `dψp(ψ)/dψ`.
    ///
    /// It's a good check that the values coincide with `qfactor.iota_of_psi(psi)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dpsip_dpsi = qfactor.dpsip_dpsi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the derivative `dψ(ψp)/dψp`.
    ///
    /// It's a good check that the values coincide with `qfactor.q_of_psip(psip)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dpsi_dpsip = qfactor.dpsi_dpsip(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `i(ψ) = 1/q(ψ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let iota_of_psi = qfactor.iota_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation of `self.q_of_psi()` fails for any reason.
    #[inline]
    fn iota_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        Ok(self.q_of_psi(psi, acc)?.recip())
    }

    /// Calculates `i(ψp) = 1/q(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let iota_of_psip = qfactor.iota_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation of `self.q_of_psip()` fails for any reason.
    #[inline]
    fn iota_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        Ok(self.q_of_psip(psip, acc)?.recip())
    }

    /// Calculates `ψ(q)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psi_of_q = qfactor.psi_of_q(1.2, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn psi_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `ψp(q)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psip_of_q = qfactor.psip_of_q(1.2, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn psip_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;
}

/// Plasma current related quantities computation.
pub trait Current {
    /// Returns the [`FluxCoordinateState`] of the toroidal `ψ` flux coordinate.
    fn psi_state(&self) -> FluxCoordinateState;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψp` flux coordinate.
    fn psip_state(&self) -> FluxCoordinateState;

    /// Calculates `g(ψ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let g_of_psi = current.g_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn g_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `g(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let g_of_psip = current.g_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn g_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `I(ψ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let i_of_psi = current.i_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn i_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `I(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let i_of_psip = current.i_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn i_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `𝜕g(ψ)/𝜕ψ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dg_dpsi = current.dg_dpsi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dg_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `𝜕g(ψp)/𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dg_dpsip = current.dg_dpsip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `𝜕I(ψ)/𝜕ψ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let di_dpsi = current.di_dpsi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn di_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates `𝜕I(ψp)/𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let di_dpsip = current.di_dpsip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError>;
}

/// Magnetic field related quantities computation.
pub trait Bfield {
    /// Returns the [`FluxCoordinateState`] of the toroidal `ψ` flux coordinate.
    fn psi_state(&self) -> FluxCoordinateState;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψp` flux coordinate.
    fn psip_state(&self) -> FluxCoordinateState;

    /// Calculates `B(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let b_of_psi = bfield.b_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn b_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `B(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let b_of_psip = bfield.b_of_psip(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn b_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `dB(ψ, θ)/dψ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_dpsi = bfield.db_dpsi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn db_dpsi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `dB(ψp, θ)/dψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_dpsip = bfield.db_dpsip(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn db_dpsip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `dB(ψ, θ)/dθ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_of_psi_dtheta = bfield.db_of_psi_dtheta(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn db_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;

    /// Calculates `dB(ψp, θ)/dθ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_of_psip_dtheta = bfield.db_of_psip_dtheta(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn db_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError>;
}

/// Defines the behavior of objects that support caching of a [`Mode`]'s values.
///
/// All non-trivial [`Mode`] evaluation methods must be guarded by a:
/// ```notest
/// if !cache.is_updated() {
///     cache.update()
/// }
/// ```
/// statement before using the cache.
#[expect(
    private_bounds,
    reason = "only used internally for creating Perturbation"
)]
pub trait ModeCache: DynModeCacheClone {
    /// Checks if the cache's stored independent coordinates are up-to-date, i.e. are equal to the
    /// passed arguments.
    fn is_updated(&mut self, flux: f64, theta: f64, zeta: f64, t: f64) -> bool;

    /// Updates the cache's coordinates and intermediate values.
    fn update(&mut self, flux: f64, theta: f64, zeta: f64, t: f64);

    /// Returns a mutable reference to the cache's `cache` array.
    fn cache(&mut self) -> &mut [f64];

    /// Returns a reference to the cache's `params` array.
    fn params(&mut self) -> &[f64];

    /// Returns a mutable reference to the cache's [`Accelerator`], if it exists.
    fn acc(&mut self) -> Option<&mut Accelerator>;

    /// Returns the cache's hits.
    fn hits(&self) -> usize;

    /// Returns the cache's misses.
    fn misses(&self) -> usize;
}

/// Single perturbation mode related quantities computation.
#[expect(
    private_bounds,
    reason = "only used internally for creating Perturbation"
)]
pub trait Mode: DynModeClone {
    /// Returns the [`FluxCoordinateState`] of the toroidal `ψ` flux coordinate.
    fn psi_state(&self) -> FluxCoordinateState;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψp` flux coordinate.
    fn psip_state(&self) -> FluxCoordinateState;

    /// Returns a default instance of the Mode's corresponding caching object.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let cache1 = mode.generate_cache();
    /// # Ok::<_, EqError>(())
    /// ```
    fn generate_cache(&self) -> DynModeCache;

    /// Calculates the mode's amplitude `α(ψ, θ, ζ, t)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let a = mode.alpha_of_psi(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn alpha_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's amplitude `α(ψp, θ, ζ, t)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let a = mode.alpha_of_psip(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn alpha_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's phase `φ(ψ, θ, ζ, t)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let phase = mode.phase_of_psi(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn phase_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's phase `φ(ψ, θ, ζ, t)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let phase = mode.phase_of_psi(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn phase_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the full mode's value `h(ψ, θ, ζ, t)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let h = mode.h_of_psi(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn h_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the full mode's value `h(ψp, θ, ζ, t)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let h = mode.h_of_psip(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn h_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to ψ, `dh(ψ, θ, ζ, t)/dψ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_dpsi = mode.dh_dpsi(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_dpsi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to ψp, `dh(ψp, θ, ζ, t)/dψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_dpsip = mode.dh_dpsip(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to θ, `dh(ψ, θ, ζ, t)/dθ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_of_psi_dtheta = mode.dh_of_psi_dtheta(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to θ, `dh(ψp, θ, ζ, t)/dθ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_of_psip_dtheta = mode.dh_of_psip_dtheta(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to ζ, `dh(ψ, θ, ζ, t)/dζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_of_psi_dzeta = mode.dh_of_psi_dzeta(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_of_psi_dzeta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to ζ, `dh(ψp, θ, ζ, t)/dζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_of_psip_dzeta = mode.dh_of_psip_dzeta(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_of_psip_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to t, `dh(ψ, θ, ζ, t)/dt`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_of_psi_dt = mode.dh_of_psi_dt(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_of_psi_dt(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to t, `dh(ψp, θ, ζ, t)/dt`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let mut cache = mode.generate_cache();
    /// let dh_of_psip_dt = mode.dh_of_psip_dt(0.1, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    fn dh_of_psip_dt(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;
}

// ================================================================================================

// HACK: This is necessary to clone `dynMode` = `Box<dyn Mode>`.
// https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object

/// Clone `Box<dyn Mode>`.
trait DynModeClone {
    /// Clone box.
    fn clone_box(&self) -> DynMode;
}

/// Clone `Box<dyn ModeCache>`.
trait DynModeCacheClone {
    /// Clone box.
    fn clone_box(&self) -> DynModeCache;
}

impl<M> DynModeClone for M
where
    M: 'static + Mode + Clone,
{
    fn clone_box(&self) -> DynMode {
        Box::new(self.clone())
    }
}

impl<C> DynModeCacheClone for C
where
    C: 'static + ModeCache + Clone,
{
    fn clone_box(&self) -> DynModeCache {
        Box::new(self.clone())
    }
}

impl Clone for DynMode {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

impl Clone for DynModeCache {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}
