//! Representation of a numerical flute mode.

use ndarray::Array1;
use rsl_interpolation::Accelerator;
use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::constants::NC_ANALYTICAL_THRESHOLD_INDEX;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{DynModeCache, Interpolation1dType, MachineObject, MachineType, Mode, ModeCache};
use crate::{EvalError, MachineError};
use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    interp_type_getter_impl, mode_cache_getters_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl,
};
use dexter_common::{DynInterpolator, make_interp};

/// Defines the calculation method of the phase `φ` in an [`NcFluteMode`].
#[derive(Default, Debug, Clone)]
pub enum PhaseMethod {
    /// Corresponds to `φ = 0`.
    #[default]
    Zero,
    /// Corresponds to `φ = const = the average of all the values of the phase array`.
    Average,
    /// Interpolation over the phase array.
    Interpolation,
    /// Use a custom value for `φ = const`.
    Custom(f64),
}

/// Used to create an [`NcFluteMode`].
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct NcFluteModeBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// The interpolation type.
    interp_type: Interpolation1dType,
    /// The `θ` frequency number.
    m: i64,
    /// The `θ` frequency number.
    n: i64,
    /// The calculation method of the phase `φ`.
    phase_method: PhaseMethod,
    /// The index of the data point under witch to use the analytical form.
    analytical_threshold_index: usize,
}

impl NcFluteModeBuilder {
    /// Creates a new [`NcFluteModeBuilder`] from a netCDF file at `path`, with spline of
    /// `interp1d_type` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcFluteModeBuilder::new(&path, Interpolation1dType::Cubic, 3, 2);
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: Interpolation1dType, m: i64, n: i64) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type,
            m,
            n,
            phase_method: PhaseMethod::default(),
            analytical_threshold_index: NC_ANALYTICAL_THRESHOLD_INDEX,
        }
    }

    /// Sets the phase `φ` calculation method.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcFluteModeBuilder::new(&path, Interpolation1dType::Cubic, 3, 2)
    ///     .with_phase_method(PhaseMethod::Interpolation)
    ///     .build()?;
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn with_phase_method(mut self, method: PhaseMethod) -> Self {
        self.phase_method = method;
        self
    }

    /// Sets the modes's analytical threshold point.
    ///
    /// By definition, flute modes must behave like `sqrt(ψ)` close to the axis, and therefore their
    /// derivative with respect to the flux must go to infinity. This is a behavior that splines
    /// cannot replicate, resulting to unnatural orbits close to the magnetic axis.
    ///
    /// To solve this, the mode switches to an analytical formula for the values of `ψ/ψp` under
    /// a certain threshold. The threshold is defined by the flux value at the position `index` of
    /// the data array.
    ///
    /// # Formula
    ///
    /// The patch has the form `β*sqrt(ψ) + γ`, where `β` and `γ` are adjusted in order to ensure
    /// continuity of both `α(ψ)` and its first derivative. `β` is calculated first by `β = 2a' * sqrt(ψ)`
    /// to ensure the correct value of the derivative `α'` at the patch's edge. Finally,
    /// `γ = α - β * sqrt(ψ)` ensures the continuity of `α` itself.
    ///
    /// Note that sometimes `γ` may become slightly negative, resulting to `α` becoming slightly
    /// negative extremely close to the axis. However this error should be negligible compared to
    /// the possible non-continuity of `α`'s higher derivatives and/or its deviation from the
    /// actual data.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let typ = Interpolation1dType::Akima;
    /// let builder = NcFluteModeBuilder::new(&path, typ, 3, 2)
    ///     .with_analytical_threshold_index(4)
    ///     .build()?;
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn with_analytical_threshold_index(mut self, index: usize) -> Self {
        self.analytical_threshold_index = index;
        self
    }

    /// Creates a new [`NcFluteMode`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let typ = Interpolation1dType::Akima;
    /// let mode = NcFluteModeBuilder::new(&path, typ, 3, 2).build()?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`MachineError`] if it fails to build the [`NcFluteMode`].
    pub fn build(self) -> Result<NcFluteMode, MachineError> {
        NcFluteMode::build(self)
    }
}

/// Single perturbation flute mode from a netCDF file.
///
/// The mode has the form of `α(ψ/ψp) * cos(mθ-nζ+φ(ψ/ψp))`, where `α` and `φ` can be expressed
/// as functions of either or both `ψ`, `ψp`, and are calculated by interpolation over the
/// numerical data.
///
/// `φ` calculation can be further configured with the [`PhaseMethod`] helper struct.
///
/// Should be created with an [`NcFluteModeBuilder`].
#[non_exhaustive]
#[derive(Clone)]
pub struct NcFluteMode {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,
    /// The interpolation type.
    interp_type: Interpolation1dType,

    /// The modes's poloidal mode number `m`.
    m: i64,
    /// The modes's toroidal mode number `n`.
    n: i64,

    /// The mode as a function of `ψ`.
    psi_single: SingleNcFluteMode,
    /// The mode as a function of `ψp`.
    psip_single: SingleNcFluteMode,
}

impl NcFluteMode {
    /// Constructs an `NcFluteMode` from an [`NcFluteModeBuilder`].
    #[expect(clippy::needless_pass_by_value, reason = "should be consumed")]
    pub(crate) fn build(builder: NcFluteModeBuilder) -> Result<Self, MachineError> {
        use crate::extract;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path.clone())?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);

        let psi_single = SingleNcFluteMode::build(&file, &builder, psi)?;
        let psip_single = SingleNcFluteMode::build(&file, &builder, psip)?;

        Ok(Self {
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            m: builder.m,
            n: builder.n,
            psi_single,
            psip_single,
        })
    }

    /// Checks if the called evaluation method is defined, returning `Err()` if not.
    fn check_if_defined(state: &FluxCoordinateState, msg: &str) -> Result<(), EvalError> {
        if *state == FluxCoordinateState::Good {
            Ok(())
        } else {
            Err(EvalError::UndefinedEvaluation(msg.into()))
        }
    }
}

impl MachineObject for NcFluteMode {
    fn machine_type(&self) -> MachineType {
        MachineType::Numerical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        self.psi_single.flux.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip_single.flux.state()
    }
}

impl std::fmt::Debug for NcFluteMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcFluteMode")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("interpolation type", &self.interp_type())
            .field("m", &self.m)
            .field("n", &self.n)
            .field("psi single", &self.psi_single)
            .field("psip single", &self.psip_single)
            .field("phase method", &self.psi_single.phase_method)
            .field(
                "analytical threshold index",
                &self.psi_single.analytical_threshold_index,
            )
            .finish()
    }
}

/// Stores an [`NcFluteMode`]'s constant parameters and cached quantities.
#[derive(Debug, Clone)]
pub struct NcFluteModeCache {
    /// The number of cache hits.
    hits: usize,
    /// The number of cache misses.
    misses: usize,
    /// Ordered array with the mode's static parameters.
    ///
    /// params = [
    ///     0 = m
    ///     1 = n
    /// ].
    params: [f64; 2],
    /// Ordered array with the modes intermediate cached values.
    ///
    /// cache = [
    ///     0 = flux
    ///     1 = theta
    ///     2 = zeta
    ///     3 = alpha
    ///     4 = dalpha
    ///     5 = phase
    ///     6 = dphase
    ///     7 = modarg
    ///     8 = sin
    ///     9 = cos
    /// ].
    cache: [f64; 10],
    /// The Accelerator of the current flux coordinate.
    acc: Accelerator,
}

impl Default for NcFluteModeCache {
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            params: [f64::NAN; 2],
            cache: [f64::NAN; 10],
            acc: Accelerator::new(),
        }
    }
}

impl ModeCache for NcFluteModeCache {
    fn is_updated(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) -> bool {
        #[expect(
            clippy::float_cmp,
            reason = "we want a cache hit only if all values are exactly equal"
        )]
        if (self.cache[0] == flux) && (self.cache[1] == theta) && (self.cache[2] == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) {
        self.cache[0] = flux;
        self.cache[1] = theta;
        self.cache[2] = zeta;

        self.cache[7] =
            (self.params()[0] * theta - self.params()[1] * zeta + self.cache[5]).rem_euclid(TAU);
        (self.cache[8], self.cache[9]) = self.cache[7].sin_cos();
    }

    fn acc(&mut self) -> Option<&mut Accelerator> {
        Some(&mut self.acc)
    }

    mode_cache_getters_impl!(NcFluteModeCache);
}

// Perform psi/psip debug assertions here since he have the extra information about the flux.
impl Mode for NcFluteMode {
    fn psi_last(&self) -> Option<f64> {
        self.psi_single.flux.last_value()
    }

    fn psip_last(&self) -> Option<f64> {
        self.psip_single.flux.last_value()
    }

    fn m(&self) -> i64 {
        self.m
    }

    fn n(&self) -> i64 {
        self.n
    }

    fn generate_cache(&self) -> DynModeCache {
        Box::new(NcFluteModeCache {
            params: [self.m as f64, self.n as f64],
            cache: [f64::NAN; 10],
            ..Default::default()
        })
    }

    fn ampl_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "α(ψ)")?;
        self.psi_single.alpha(psi, theta, zeta, t, cache)
    }

    fn ampl_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "α(ψp)")?;
        self.psip_single.alpha(psip, theta, zeta, t, cache)
    }

    fn phase_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "φ(ψ)")?;
        self.psi_single.phase(psi, theta, zeta, t, cache)
    }

    fn phase_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "φ(ψp)")?;
        self.psip_single.phase(psip, theta, zeta, t, cache)
    }

    fn m_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "h(ψ)")?;
        self.psi_single.m(psi, theta, zeta, t, cache)
    }

    fn m_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "h(ψp)")?;
        self.psip_single.m(psip, theta, zeta, t, cache)
    }

    fn dm_dpsi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "dh(ψ)/dψ")?;
        self.psi_single.dm_dflux(psi, theta, zeta, t, cache)
    }

    fn dm_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "dh(ψp)/dψp")?;
        self.psip_single.dm_dflux(psip, theta, zeta, t, cache)
    }

    fn dm_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "dh(ψ)/dθ")?;
        self.psi_single.dm_dtheta(psi, theta, zeta, t, cache)
    }

    fn dm_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "dh(ψp)/dθ")?;
        self.psip_single.dm_dtheta(psip, theta, zeta, t, cache)
    }

    fn dm_of_psi_dzeta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "dh(ψ)/dζ")?;
        self.psi_single.dm_dzeta(psi, theta, zeta, t, cache)
    }

    fn dm_of_psip_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "dh(ψp)/dζ")?;
        self.psip_single.dm_dzeta(psip, theta, zeta, t, cache)
    }

    fn dm_of_psi_dt(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(0.0)
    }

    fn dm_of_psip_dt(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(0.0)
    }
}

// Getters.
impl NcFluteMode {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    interp_type_getter_impl!(1);

    /// Returns the mode's [`PhaseMethod`].
    ///
    /// Both flux coordinates have the same [`PhaseMethod`].
    #[must_use]
    pub fn phase_method(&self) -> PhaseMethod {
        self.psi_single.phase_method.clone()
    }

    /// Returns the average value of the phase arrays, if [`PhaseMethod`] is `Average`.
    #[must_use]
    pub fn phase_average(&self) -> Option<f64> {
        self.psi_single.phase_average
    }

    /// Returns the toroidal flux's values as a 1D array, if they exist.
    #[must_use]
    pub fn psi_array(&self) -> Option<Array1<f64>> {
        self.psi_single
            .flux
            .values()
            .map(|values| Array1::from(Vec::from(values)))
    }
    /// Returns the poloidal flux's values as a 1D array, if they exist.
    #[must_use]
    pub fn psip_array(&self) -> Option<Array1<f64>> {
        self.psip_single
            .flux
            .values()
            .map(|values| Array1::from(Vec::from(values)))
    }

    /// Returns the `α` values as a 1D array.
    #[must_use]
    pub fn alpha_array(&self) -> Array1<f64> {
        Array1::from(self.psi_single.alpha_values.clone())
    }

    /// Returns the `φ` values as a 1D array.
    #[must_use]
    pub fn phase_array(&self) -> Array1<f64> {
        Array1::from(self.psi_single.phase_values.clone())
    }

    /// Returns the index of the analytical threshold.
    #[must_use]
    pub fn analytical_threshold_index(&self) -> usize {
        // This exists in both single modes
        self.psi_single.analytical_threshold_index
    }
}

// ===============================================================================================
// ===============================================================================================

/// Representation of a numerical flute mode, defined as a function of only one of the two flux
/// coordinates (and θ, ζ, t).
///
/// Splitting an [`NcFluteMode`] into two `SingleNcFluteModes` makes the code much clearer.
///
/// Both modes end up identical, with the only difference being the corresponding [`NcFlux`].
#[non_exhaustive]
#[derive(Clone)]
struct SingleNcFluteMode {
    /// The current flux coordinate.
    flux: NcFlux,
    /// The phase calculation method.
    phase_method: PhaseMethod,
    /// The phase values' average, if `phase_method` is `Average`.
    phase_average: Option<f64>,
    /// The index of the data point under witch to use the analytical formula.
    analytical_threshold_index: usize,
    /// The flux's value at the `analytical_threshold_index`.
    analytical_threshold_flux: Option<f64>,
    /// The patch's `β` coefficient.
    patch_beta: Option<f64>,
    /// The patch's `γ` coefficient.
    patch_gamma: Option<f64>,

    /// The amplitude values.
    alpha_values: Box<[f64]>,
    /// The phase values.
    phase_values: Box<[f64]>,
    /// The `α(flux)` interpolator.
    alpha_interp: Option<DynInterpolator>,
    /// The `φ(flux)` interpolator.
    phase_interp: Option<DynInterpolator>,
}

// Creation.
impl SingleNcFluteMode {
    /// Creates a `SingleNcFluteMode` from a [`NcFlux`].
    ///
    /// The object always exists, even if the corresponding flux coordinate does not exist, so all
    /// logic and error handling is handle here.
    pub(crate) fn build(
        file: &netcdf::File,
        builder: &NcFluteModeBuilder,
        flux: NcFlux,
    ) -> Result<Self, MachineError> {
        use crate::extract;

        let (alpha_data, phase_data) = extract::mode_arrays(file, builder.m, builder.n)?;
        let alpha_values = alpha_data.to_vec().into_boxed_slice();
        let phase_values = phase_data.to_vec().into_boxed_slice();

        debug_assert_all_finite_values(&alpha_values);
        debug_assert_all_finite_values(&phase_values);

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let alpha_interp = match flux.state() {
            Good => Some(make_interp(
                builder.interp_type,
                flux.uvalues(),
                &alpha_values,
            )?),
            _ => None,
        };
        #[rustfmt::skip]
        let phase_interp = match flux.state() {
            Good => Some(make_interp(
                builder.interp_type,
                flux.uvalues(),
                &phase_values,
            )?),
            _ => None,
        };

        let _m = builder.m as f64;
        let _n = builder.n as f64;

        let self_state0 = Self {
            flux,
            phase_method: builder.phase_method.clone(),
            phase_average: None,
            analytical_threshold_index: builder.analytical_threshold_index,
            analytical_threshold_flux: None,
            patch_beta: None,
            patch_gamma: None,
            alpha_values,
            phase_values,
            alpha_interp,
            phase_interp,
        };
        let self_state1 = self_state0.resolve_phase_method();
        let final_state = self_state1.patch_axis()?;
        Ok(final_state)
    }

    /// Calculates the phase method. Should only be used on initialization.
    fn resolve_phase_method(self) -> Self {
        let phase_method: PhaseMethod;
        let mut phase_average: Option<f64> = None;

        match self.phase_method {
            PhaseMethod::Zero => phase_method = PhaseMethod::Zero,
            PhaseMethod::Interpolation => phase_method = PhaseMethod::Interpolation,
            PhaseMethod::Custom(phase) => phase_method = PhaseMethod::Custom(phase),
            PhaseMethod::Average => {
                phase_method = PhaseMethod::Average;
                phase_average = Array1::from(self.phase_values.clone()).mean();
            }
        }

        Self {
            phase_method,
            phase_average,
            ..self
        }
    }

    /// Calculates the patch's `β` and `γ` coefficients.
    #[expect(clippy::unwrap_in_result, reason = "cannot panic")]
    fn patch_axis(mut self) -> Result<Self, MachineError> {
        let Some(interp) = self.alpha_interp.as_ref() else {
            return Ok(self);
        };

        let acc = &mut Accelerator::new();
        let flux_value = self
            .flux
            .uvalues()
            .get(self.analytical_threshold_index)
            .copied()
            .ok_or(MachineError::InvalidFluteModeAnalyticalThresholdIndex)?;
        let switch_alpha = interp
            .eval(self.flux.uvalues(), &self.alpha_values, flux_value, acc)
            .expect("domain just checked");
        let switch_dalpha = interp
            .eval_deriv(self.flux.uvalues(), &self.alpha_values, flux_value, acc)
            .expect("domain just checked");

        // 1st derivative continuity condition
        let patch_beta = 2.0 * switch_dalpha * flux_value.sqrt();
        // `α` continuity condition
        let patch_gamma = switch_alpha - patch_beta * flux_value.sqrt();

        self.analytical_threshold_flux = Some(flux_value);
        self.patch_beta = Some(patch_beta);
        self.patch_gamma = Some(patch_gamma);

        Ok(self)
    }
}

/// External cache update.
impl SingleNcFluteMode {
    /// Updates the interpolated values, since they cannot take place inside the cache without
    /// messing up the whole structure.
    ///
    /// Should always be called together with `cache.update()`, and `update_cache_interps()` should
    /// always be called firsts, since `cache` uses the phase value.
    #[rustfmt::skip]
    fn update_cache_interps(
        &self,
        flux: f64,
        cache: &mut DynModeCache,
    ) -> Result<(), EvalError> {

        if self.analytical_threshold_flux.is_some_and(|threshold| flux < threshold) {
            let root = flux.sqrt();
            let patch_beta = self.patch_beta.expect("just checked");
            let patch_gamma = self.patch_gamma.expect("just checked");
            cache.cache()[3] = patch_beta * root + patch_gamma;
            cache.cache()[4] = patch_beta / (2.0 * root);
        } else {
            cache.cache()[3] = match self.alpha_interp.as_ref() {
                Some(interp) => interp.eval(self.flux.uvalues(), &self.alpha_values, flux, cache.acc().expect("has"))?,
                None => return Err(EvalError::UndefinedEvaluation("α(flux)".into())),
            };
            cache.cache()[4] = match self.alpha_interp.as_ref() {
                Some(interp) => interp.eval_deriv(self.flux.uvalues(), &self.alpha_values, flux, cache.acc().expect("has"))?,
                None => return Err(EvalError::UndefinedEvaluation("da(flux)/dflux".into())),
            };
        }

        cache.cache()[5] = match self.phase_method {
            PhaseMethod::Zero => 0.0,
            PhaseMethod::Average => self.phase_average.expect("Exists"),
            PhaseMethod::Custom(custom_phase) => custom_phase,
            PhaseMethod::Interpolation => {
                match self.phase_interp.as_ref() {
                    Some(interp) => interp.eval(self.flux.uvalues(), &self.phase_values, flux, cache.acc().expect("has"))?,
                    None => return Err(EvalError::UndefinedEvaluation("φ(flux)".into())),
                }
            }
        };
        cache.cache()[6] = match self.phase_method {
            PhaseMethod::Interpolation => {
                match self.phase_interp.as_ref() {
                    Some(interp) => interp.eval_deriv(self.flux.uvalues(), &self.phase_values, flux, cache.acc().expect("has"))?,
                    None => return Err(EvalError::UndefinedEvaluation("φ(flux)".into())),
                }
            }
            _ => 0.0,
        };
        Ok(())
    }
}

/// Intermediate Interpolations. This is effectively where the [`Mode`] trait is implemented.
impl SingleNcFluteMode {
    /// Calculates the single modes's amplitude.
    fn alpha(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache()[3]))
    }

    /// Calculates the single mode's phase.
    fn phase(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache()[5]))
    }

    /// Calculates the single mode's value.
    fn m(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache()[3] * cache.cache()[9]))
    }

    /// Calculates the single mode's derivative with respect to the current flux coordinate.
    fn dm_dflux(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache()[4] * cache.cache()[9])
            - cache.cache()[3] * cache.cache()[8] * cache.cache()[6])
    }

    /// Calculates the single mode's derivative with respect to theta.
    fn dm_dtheta(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(
            -cache.params()[0] * cache.cache()[3] * cache.cache()[8]
        ))
    }

    /// Calculates the single mode's derivative with respect to zeta.
    fn dm_dzeta(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(
            cache.params()[1] * cache.cache()[3] * cache.cache()[8]
        ))
    }
}

impl std::fmt::Debug for SingleNcFluteMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SingleNcFluteMode")
            .field("NcFlux", &self.flux)
            .field("phase method", &self.phase_method)
            .field("phase average", &self.phase_average)
            .field(
                "analytical threshold",
                &format!("{:?}", self.analytical_threshold_flux),
            )
            .field(
                "analytical threshold patch-beta",
                &format!("{:?}", self.patch_beta),
            )
            .field(
                "analytical threshold patch-gamma",
                &format!("{:?}", self.patch_gamma),
            )
            .finish()
    }
}

// ===============================================================================================

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_flute_mode_builder(path: &str) -> NcFluteModeBuilder {
        let path = PathBuf::from(path);
        NcFluteModeBuilder::new(&path, Interpolation1dType::Cubic, 3, 2)
    }

    pub(super) fn create_nc_flute_mode(path: &str) -> NcFluteMode {
        create_nc_flute_mode_builder(path)
            .with_phase_method(PhaseMethod::Zero)
            .build()
            .unwrap()
    }
}

#[cfg(test)]
mod phase_methods {
    use crate::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
    use approx::assert_relative_eq;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn zero_phase_method() {
        use PhaseMethod::Zero;
        let mode = create_nc_flute_mode_builder(TEST_NETCDF_PATH)
            .with_phase_method(Zero)
            .build()
            .unwrap();
        let c = &mut mode.generate_cache();

        assert!(matches!(mode.phase_method(), Zero));
        assert!(mode.psi_single.phase_average.is_none());

        assert_eq!(mode.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), 0.0);
        assert_eq!(mode.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), 0.0);
    }

    #[test]
    fn average_phase_method() {
        use PhaseMethod::Average;
        let mode = create_nc_flute_mode_builder(TEST_NETCDF_PATH)
            .with_phase_method(Average)
            .build()
            .unwrap();
        let c = &mut mode.generate_cache();

        let expected = mode.phase_array().mean().unwrap();

        assert!(matches!(mode.phase_method(), Average));
        assert!(mode.psi_single.phase_average.is_some());

        assert_eq!(mode.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), expected);
        assert_eq!(
            mode.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(),
            expected
        );
    }

    #[test]
    fn interpolation_phase_method() {
        use PhaseMethod::Interpolation;
        let mode = create_nc_flute_mode_builder(POLOIDAL_TEST_NETCDF_PATH)
            .with_phase_method(Interpolation)
            .build()
            .unwrap();
        let c = &mut mode.generate_cache();

        // Calculated with a stable version, on the same dataset
        let expected = 0.8414720460746888;

        assert!(matches!(mode.phase_method(), Interpolation));
        assert!(mode.psi_single.phase_average.is_none());

        assert_relative_eq!(
            mode.phase_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap(),
            expected,
            epsilon = 1e-5 // unsure why there is this small difference here.
        );
    }

    #[test]
    fn custom_phase_method() {
        use PhaseMethod::Custom;
        let mode = create_nc_flute_mode_builder(TEST_NETCDF_PATH)
            .with_phase_method(Custom(10.0))
            .build()
            .unwrap();
        let c = &mut mode.generate_cache();

        assert!(matches!(mode.phase_method(), Custom(10.0)));
        assert!(mode.psi_single.phase_average.is_none());

        assert_eq!(mode.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), 10.0);
        assert_eq!(mode.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), 10.0);
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let mode = create_nc_flute_mode(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(mode.psi_state(), FluxCoordinateState::Good);
        assert_eq!(mode.psip_state(), FluxCoordinateState::Bad);
        assert!(mode.psi_single.alpha_interp.is_some());
        assert!(mode.psi_single.phase_interp.is_some());
        assert!(mode.psip_single.alpha_interp.is_none());
        assert!(mode.psip_single.phase_interp.is_none());

        assert!(mode.psi_array().is_some());
        assert!(mode.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psi_evals() {
        let mode = create_nc_flute_mode(TOROIDAL_TEST_NETCDF_PATH);
        let c = &mut mode.generate_cache();
        assert!(mode.ampl_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.phase_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.m_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_dpsi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_of_psi_dtheta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_of_psi_dzeta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_of_psi_dt(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
    }

    #[test]
    #[rustfmt::skip]
    fn bad_psip_evals() {
        let mode = create_nc_flute_mode(TOROIDAL_TEST_NETCDF_PATH);
        let c = &mut mode.generate_cache();

        use EvalError::UndefinedEvaluation as err;
        assert!(matches!(mode.ampl_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.phase_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.m_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_dpsip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_of_psip_dtheta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_of_psip_dzeta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_of_psip_dt(0.1, 0.1, 0.1, 0.1, c), Ok(0.0)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let mode = create_nc_flute_mode(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(mode.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(mode.psip_state(), FluxCoordinateState::Good);
        assert!(mode.psi_single.alpha_interp.is_none());
        assert!(mode.psi_single.phase_interp.is_none());
        assert!(mode.psip_single.alpha_interp.is_some());
        assert!(mode.psip_single.phase_interp.is_some());

        assert!(mode.psi_array().is_some());
        assert!(mode.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psip_evals() {
        let mode = create_nc_flute_mode(POLOIDAL_TEST_NETCDF_PATH);
        let c = &mut mode.generate_cache();
        assert!(mode.ampl_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.phase_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.m_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_dpsip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_of_psip_dtheta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_of_psip_dzeta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(mode.dm_of_psip_dt(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
    }

    #[test]
    #[rustfmt::skip]
    fn bad_psi_evals() {
        let mode = create_nc_flute_mode(POLOIDAL_TEST_NETCDF_PATH);
        let c = &mut mode.generate_cache();

        use EvalError::UndefinedEvaluation as err;
        assert!(matches!(mode.ampl_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.phase_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.m_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_dpsi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_of_psi_dtheta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_of_psi_dzeta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(mode.dm_of_psi_dt(0.1, 0.1, 0.1, 0.1, c), Ok(0.0)));
    }
}

#[cfg(test)]
mod nc_flute_mode_cache {
    use crate::extract::TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    #[allow(unused_results)]
    fn counts() {
        let mode = create_nc_flute_mode_builder(TEST_NETCDF_PATH)
            // other methods don't touch the cache at all
            .with_phase_method(PhaseMethod::Interpolation)
            .build()
            .unwrap();
        let c = &mut mode.generate_cache();
        let t = 0.0; // not checked

        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 0);

        mode.phase_of_psi(0.01, 0.1, 0.1, t, c).unwrap();
        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 1);

        mode.dm_of_psi_dtheta(0.01, 0.1, 0.1, t, c).unwrap();
        assert_eq!(c.hits(), 1);
        assert_eq!(c.misses(), 1);

        let psi = 0.1;
        let theta = 3.14;
        let zeta = 1.0;
        mode.ampl_of_psi(psi, theta, zeta, t, c).unwrap();
        mode.phase_of_psi(psi, theta, zeta, t, c).unwrap();
        mode.m_of_psi(psi, theta, zeta, t, c).unwrap();
        mode.dm_of_psi_dtheta(psi, theta, zeta, t, c).unwrap();
        mode.dm_of_psi_dzeta(psi, theta, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 5);
        assert_eq!(c.misses(), 2);

        mode.dm_dpsi(psi / 2.0, theta, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 5);
        assert_eq!(c.misses(), 3);

        // dbg!(&c); // Accelerator counts should also match
    }
}

#[cfg(test)]
mod nc_flute_mode_analytical_threshold {
    use crate::extract::TEST_NETCDF_PATH;

    use super::*;

    #[test]
    fn normal_construction() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let mode = NcFluteModeBuilder::new(&path, Interpolation1dType::Cubic, 3, 2)
            .with_phase_method(PhaseMethod::Zero)
            .with_analytical_threshold_index(5)
            .build()
            .unwrap();

        let mut cache = mode.generate_cache();
        assert_eq!(mode.analytical_threshold_index(), 5);
        // set cos=1 and make sure it does go to infinity
        assert!(
            dbg!(mode.dm_dpsi(1e-20, 0.0, 0.0, 0.0, &mut cache))
                .is_ok_and(|value| value.abs() > 1000.0)
        )
    }

    #[test]
    fn construction_with_zero_index() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let mode_attempt = NcFluteModeBuilder::new(&path, Interpolation1dType::Cubic, 3, 2)
            .with_phase_method(PhaseMethod::Zero)
            .with_analytical_threshold_index(0)
            .build();
        assert!(mode_attempt.is_ok());

        let mode = mode_attempt.unwrap();

        let mut cache = mode.generate_cache();
        assert!(
            mode.m_of_psi(1e-5, 0.0, 0.0, 0.0, &mut cache)
                .unwrap()
                .is_finite()
        );
    }

    #[test]
    fn erroneous_construction() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let mode_attempt = NcFluteModeBuilder::new(&path, Interpolation1dType::Cubic, 3, 2)
            .with_phase_method(PhaseMethod::Zero)
            .with_analytical_threshold_index(5000000000)
            .build();

        assert!(mode_attempt.is_err_and(|err| matches!(
            err,
            MachineError::InvalidFluteModeAnalyticalThresholdIndex
        )));
    }
}
