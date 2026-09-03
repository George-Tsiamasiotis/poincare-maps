//! Representation of a machines's q-factor profile.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    fluxes_values_array_getter_impl, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl,
};
use ndarray::Array1;
use rsl_interpolation::Accelerator;
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::Interpolation1dType;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{EvalError, MachineError};
use crate::{FluxCommute, LastClosedFluxSurface, MachineObject, MachineType, Qfactor};
use dexter_common::{DynInterpolator, array1D_getter_impl, make_interp};

// ===============================================================================================

/// Analytical q-factor profile of q = 1 and ψ=ψp.
#[non_exhaustive]
#[derive(Clone)]
pub struct UnityQfactor {
    /// The value of the last closed toroidal flux surface `ψ_last` in Normalized units.
    psi_last: f64,
    /// The value of the last closed poloidal flux surface `ψp_last` in Normalized units.
    psip_last: f64,
}

impl UnityQfactor {
    /// Creates a new `UnityQfactor`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// // Define ψ=ψ_last=0.45
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let qfactor = UnityQfactor::new(lcfs);
    /// ```
    #[must_use]
    pub fn new(lcfs: LastClosedFluxSurface) -> Self {
        Self {
            psi_last: lcfs.value(),
            psip_last: lcfs.value(),
        }
    }
}

impl MachineObject for UnityQfactor {
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

impl FluxCommute for UnityQfactor {
    fn psip_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if psi > self.psi_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(psi)
    }

    fn psi_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if psip > self.psip_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(psip)
    }
}

impl Qfactor for UnityQfactor {
    fn psi_last(&self) -> f64 {
        self.psi_last
    }

    fn psip_last(&self) -> f64 {
        self.psip_last
    }

    fn qlast(&self) -> f64 {
        1.0
    }

    fn qaxis(&self) -> f64 {
        1.0
    }

    fn q_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if psi > self.psi_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }

    fn q_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if psip > self.psip_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }

    fn dpsip_dpsi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if psi > self.psi_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }

    fn dpsi_dpsip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if psip > self.psip_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }

    fn psi_of_q(&self, _: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation("ψ(q)".into()))
    }

    fn psip_of_q(&self, _: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation("ψp(q)".into()))
    }
}

impl std::fmt::Debug for UnityQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("q-factor of q = 1 and ψ=ψp")
            .field("ψ/ψp_last", &self.psi_last)
            .finish()
    }
}

// ===============================================================================================

/// Analytical q-factor of parabolic q(ψ) profile.
#[derive(Clone)]
pub struct ParabolicQfactor {
    /// The `q` value on the magnetic axis.
    qaxis: f64,
    /// The `q` value at the last closed flux surface.
    qlast: f64,
    /// The value of the last closed toroidal flux surface `ψ_last` in Normalized units.
    psi_last: f64,
    /// The value of the last closed poloidal flux surface `ψp_last` in Normalized units.
    psip_last: f64,
}

impl ParabolicQfactor {
    /// Creates a new `ParabolicQfactor`.
    ///
    /// A `ParabolicQfactor` is defined with the help of the [`LastClosedFluxSurface`] helper struct,
    /// which changes the position where `qlast` is met.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// // Define q(ψ=ψ_last=0.45) = qlast = 3.8
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    ///
    /// // Define q(ψp=ψp_last=0.19) = qlast = 4.2
    /// let psip_last = LastClosedFluxSurface::Poloidal(0.19);
    /// let qfactor = ParabolicQfactor::new(1.1, 4.2, psip_last);
    /// ```
    #[must_use]
    pub fn new(qaxis: f64, qlast: f64, lcfs: LastClosedFluxSurface) -> Self {
        // Create two phony qfactors to calculate the other last value
        let psi_last: f64;
        let psip_last: f64;
        match lcfs {
            LastClosedFluxSurface::Toroidal(_psi_last) => {
                let phony_q = Self {
                    qaxis,
                    qlast,
                    psi_last: _psi_last,
                    psip_last: f64::NAN,
                };
                psi_last = phony_q.psi_last;
                psip_last = match phony_q.psip_of_psi(_psi_last, &mut Accelerator::new()) {
                    Ok(value) => value,
                    Err(_) => unreachable!("Analytical formula, cannot fail"),
                };
            }
            LastClosedFluxSurface::Poloidal(_psip_last) => {
                let phony_q = Self {
                    qaxis,
                    qlast,
                    psi_last: f64::NAN,
                    psip_last: _psip_last,
                };
                psip_last = phony_q.psip_last;
                psi_last = phony_q.psi_last_of_psip_last();
            }
        }

        Self {
            qaxis,
            qlast,
            psi_last,
            psip_last,
        }
    }

    /// Helper function to calculate `ψ_last` from `ψp_last` when instantiating `ParabolicQfactor`. This
    /// little maneuver is necessary since `ψ(ψp)` is defined through `ψ_last`, which of course does
    /// not exist yet.
    ///
    /// This formula is derived by solving `q(ψp)` for `ψ_last` and setting `ψp = ψp_last`.
    ///
    /// Should not be used after instantiation.
    #[must_use]
    fn psi_last_of_psip_last(&self) -> f64 {
        let atan_arg = (self.qlast / self.qaxis - 1.0).sqrt();
        let numerator = self.psip_last * (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        let denominator = atan_arg.atan();
        numerator / denominator
    }
}

impl MachineObject for ParabolicQfactor {
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

impl FluxCommute for ParabolicQfactor {
    fn psip_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if psi > self.psi_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        let atan_arg = psi * (self.qlast - self.qaxis).sqrt() / (self.psi_last * self.qaxis.sqrt());
        let coef = self.psi_last / (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        Ok(debug_assert_is_finite!(coef * atan_arg.atan()))
    }

    fn psi_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if psip > self.psip_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        let tan_arg = (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psi_last;
        let coef = self.psi_last * self.qaxis.sqrt() / (self.qlast - self.qaxis).sqrt();
        Ok(debug_assert_is_finite!(coef * tan_arg.tan()))
    }
}

impl Qfactor for ParabolicQfactor {
    fn psi_last(&self) -> f64 {
        self.psi_last
    }

    fn psip_last(&self) -> f64 {
        self.psip_last
    }

    fn qaxis(&self) -> f64 {
        self.qaxis
    }

    fn qlast(&self) -> f64 {
        self.qlast
    }

    fn q_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if psi > self.psi_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(debug_assert_is_finite!(
            (self.qaxis) + (self.qlast - self.qaxis) * (psi / self.psi_last).powi(2)
        ))
    }

    fn q_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if psip > self.psip_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        let tan_arg = (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psi_last;
        Ok(debug_assert_is_finite!(
            self.qaxis + self.qaxis * tan_arg.tan().powi(2)
        ))
    }

    fn dpsip_dpsi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if psi > self.psi_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        let denom = self.qaxis * self.psi_last.powi(2) + (self.qlast - self.qaxis) * psi.powi(2);
        Ok(debug_assert_is_finite!(self.psi_last.powi(2) / denom))
    }

    fn dpsi_dpsip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if psip > self.psip_last {
            return Err(EvalError::AnalyticalDomainError);
        };
        let cos_arg = (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psi_last;
        Ok(debug_assert_is_finite!(self.qaxis / cos_arg.cos().powi(2)))
    }

    fn psi_of_q(&self, q: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        if q > self.qlast {
            return Err(EvalError::AnalyticalDomainError);
        };
        let frac = (q - self.qaxis) / (self.qlast - self.qaxis);
        let psi = self.psi_last * frac.sqrt();
        debug_assert_non_negative_psi!(psi);
        Ok(psi)
    }

    fn psip_of_q(&self, q: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        if q > self.qlast {
            return Err(EvalError::AnalyticalDomainError);
        };
        let tan_arg = (q / self.qaxis - 1.0).sqrt();
        let coef = self.psi_last / (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        let psip = coef * tan_arg.atan();
        debug_assert_non_negative_psip!(psip);
        Ok(psip)
    }
}

impl std::fmt::Debug for ParabolicQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParabolicQfactor: q-factor of parabolic q(ψ) profile.")
            .field("ψ_last", &self.psi_last)
            .field("ψp_last", &self.psip_last)
            .field("qaxis", &self.qaxis)
            .field("qlast", &self.qlast)
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcQfactor`].
///
/// Exists for future configuration flexibility.
#[derive(Debug, Clone)]
pub struct NcQfactorBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// The interpolation type.
    interp_type: Interpolation1dType,
}

impl NcQfactorBuilder {
    /// Creates a new `NcQfactorBuilder` from a netCDF file at `path`, with `interp_type`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, Interpolation1dType::Cubic);
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: Interpolation1dType) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type,
        }
    }

    /// Creates a new [`NcQfactor`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Akima).build()?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`MachineError`] if it fails to build the [`NcQfactor`].
    pub fn build(self) -> Result<NcQfactor, MachineError> {
        NcQfactor::build(self)
    }
}

/// Numerical q-factor profile reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcQfactorBuilder`].
///
/// If either `psi_norm` or `psip_norm` is missing from the netCDF file, it is calculated from the
/// other by integrating `q(ψp)` or `ι(ψ)` respectively. In the case that the calculated values are
/// monotonic, the other flux can be used as a flux coordinate as well.
pub struct NcQfactor {
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

    /// `ψp(ψ)` interpolator.
    psip_of_psi_interp: Option<DynInterpolator>,
    /// `ψ(ψp)` interpolator.
    psi_of_psip_interp: Option<DynInterpolator>,

    /// The `q` values.
    q_values: Box<[f64]>,
    /// `q(ψ)` interpolator.
    q_of_psi_interp: Option<DynInterpolator>,
    /// `q(ψp)` interpolator.
    q_of_psip_interp: Option<DynInterpolator>,

    /// `ψp(q)` interpolator.
    psi_of_q_interp: Option<DynInterpolator>,
    /// `ψp(q)` interpolator.
    psip_of_q_interp: Option<DynInterpolator>,
}

/// Creation.
impl NcQfactor {
    /// Constructs an `NcQfactor` from an [`NcQfactorBuilder`].
    fn build(builder: NcQfactorBuilder) -> Result<Self, MachineError> {
        use crate::extract;
        use crate::extract::netcdf_fields::NC_Q;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let mut psi = NcFlux::toroidal(&file);
        let mut psip = NcFlux::poloidal(&file);
        let q_values = extract::array_1d(&file, NC_Q)?.to_vec();

        debug_assert_all_finite_values(&q_values);

        // If no `ψ` values found, integrate `q(ψp)` and create them. We create a temporary
        // `q_of_psip_interp` since we cannot correctly define it yet.
        //
        // This operation cannot fail. If one of the fluxes is missing, then the other one is
        // guaranteed to be `Good` and therefore the interpolation and integration can be defined.
        if psi.state() == FluxCoordinateState::NoValues {
            let acc = &mut Accelerator::new();
            let psip_values = psip.values().expect("At least one of the fluxes exists");
            let q_of_psip_interp = make_interp(builder.interp_type, psip_values, &q_values)?;
            let psi_values: Vec<f64> = psip_values
                .iter()
                .map(|psip_value| {
                    q_of_psip_interp
                        .eval_integ(psip_values, &q_values, 0.0, *psip_value, acc)
                        .expect("Cannot fail by definition of the interpolator")
                })
                .collect();
            psi = NcFlux::from_raw_values(&psi_values);
        }

        // If no `ψp` values found, integrate `ι(ψ)` and create them. For this we must create a
        // temporary `i_of_psi` interpolator.
        //
        // This operation cannot fail. If one of the fluxes is missing, then the other one is
        // guaranteed to be `Good` and therefore the interpolation and integration can be defined.
        if psip.state() == FluxCoordinateState::NoValues {
            let acc = &mut Accelerator::new();
            let psi_values = psi.values().expect("At least one of the fluxes exists");
            let i_values: Vec<f64> = q_values.iter().map(|q| q.recip()).collect();
            let i_of_psi_interp = make_interp(builder.interp_type, psi_values, &i_values)?;
            let psip_values: Vec<f64> = psi_values
                .iter()
                .map(|psi_value| {
                    i_of_psi_interp
                        .eval_integ(psi_values, &i_values, 0.0, *psi_value, acc)
                        .expect("Cannot fail by definition of the interpolator")
                })
                .collect();
            psip = NcFlux::from_raw_values(&psip_values);
        }

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let psip_of_psi_interp =
            if (psi.state() == Good) & (psip.state() != FluxCoordinateState::NoValues) {
                Some(make_interp(
                    builder.interp_type,
                    psi.uvalues(),
                    psip.uvalues(),
                )?)
            } else {
                None
            };
        let psi_of_psip_interp =
            if (psip.state() == Good) & (psi.state() != FluxCoordinateState::NoValues) {
                Some(make_interp(
                    builder.interp_type,
                    psip.uvalues(),
                    psi.uvalues(),
                )?)
            } else {
                None
            };

        let q_of_psi_interp = match psi.state() {
            Good => Some(make_interp(builder.interp_type, psi.uvalues(), &q_values)?),
            _ => None,
        };
        let q_of_psip_interp = match psip.state() {
            Good => Some(make_interp(builder.interp_type, psip.uvalues(), &q_values)?),
            _ => None,
        };

        // If flux values exist, we must also check if q is monotonic
        let psi_of_q_interp = match psi.state() {
            FluxCoordinateState::NoValues => None,
            _ => make_interp(builder.interp_type, &q_values, psi.uvalues()).ok(),
        };
        let psip_of_q_interp = match psip.state() {
            FluxCoordinateState::NoValues => None,
            _ => make_interp(builder.interp_type, &q_values, psip.uvalues()).ok(),
        };

        Ok(Self {
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            psi,
            psip,
            psip_of_psi_interp,
            psi_of_psip_interp,
            q_values: q_values.into_boxed_slice(),
            q_of_psi_interp,
            q_of_psip_interp,
            psi_of_q_interp,
            psip_of_q_interp,
        })
    }
}

impl MachineObject for NcQfactor {
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

impl FluxCommute for NcQfactor {
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

impl Qfactor for NcQfactor {
    fn psi_last(&self) -> f64 {
        self.psi
            .last_value()
            .expect("'psi' values have been calculated at this point")
    }

    fn psip_last(&self) -> f64 {
        self.psip
            .last_value()
            .expect("'psip' values have been calculated at this point")
    }

    fn qlast(&self) -> f64 {
        self.q_array()
            .last()
            .copied()
            .expect("'q' array is non-empty")
    }

    fn qaxis(&self) -> f64 {
        self.q_array()
            .first()
            .copied()
            .expect("'q' array is non-empty")
    }

    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.q_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.q_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("q(ψ)".into())),
        }
    }

    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.q_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.q_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("q(ψp)".into())),
        }
    }

    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.psip_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psi.uvalues(),
                self.psip.uvalues(),
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dψp(ψ)/dψ".into())),
        }
    }

    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.psi_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psip.uvalues(),
                self.psi.uvalues(),
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dψ(ψp)/dψp".into())),
        }
    }

    fn psi_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        match self.psi_of_q_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                &self.q_values,
                self.psi.uvalues(),
                q,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψ(q)".into())),
        }
    }

    fn psip_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        match self.psip_of_q_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                &self.q_values,
                self.psip.uvalues(),
                q,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψp(q)".into())),
        }
    }
}

/// Getters.
impl NcQfactor {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    interp_type_getter_impl!(1);
    fluxes_values_array_getter_impl!();
    array1D_getter_impl!(q_array, q_values, q);
}

impl std::fmt::Debug for NcQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcQfactor")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("interpolation type", &self.interp_type())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .field("qaxis", &self.qaxis())
            .field("qlast", &self.qlast())
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;
    use approx::assert_relative_eq;

    pub(super) fn create_nc_qfactor(path_str: &str) -> NcQfactor {
        let path = PathBuf::from(&path_str);
        let builder = NcQfactorBuilder::new(&path, Interpolation1dType::Steffen);
        builder.build().unwrap()
    }

    /// Make sure that dψ(ψp)/dψp and q(ψ) are close enough.
    pub(super) fn test_dpsi_dpsip_q_closeness(qfactor: &impl Qfactor) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psips = Array1::linspace(0.02 * qfactor.psip_last(), 0.98 * qfactor.psip_last(), 100);

        let mut acc = Accelerator::new();
        for psip in psips.iter().copied() {
            assert_relative_eq!(
                qfactor.q_of_psip(psip, &mut acc).unwrap(),
                qfactor.dpsi_dpsip(psip, &mut acc).unwrap(),
                epsilon = qfactor.qlast() * 1e-4
            )
        }
    }

    /// Make sure that dψp(ψ)/dψ and i(ψ) are close enough.
    pub(super) fn test_dpsip_dpsi_iota_closeness(qfactor: &impl Qfactor) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psis = Array1::linspace(0.02 * qfactor.psi_last(), 0.98 * qfactor.psi_last(), 100);

        let mut acc = Accelerator::new();
        for psi in psis.iter().copied() {
            assert_relative_eq!(
                qfactor.iota_of_psi(psi, &mut acc).unwrap(),
                qfactor.dpsip_dpsi(psi, &mut acc).unwrap(),
                epsilon = qfactor.qlast() * 1e-4
            )
        }
    }
}

#[cfg(test)]
mod test_parabolic_qfactor {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    #[test]
    /// Values calculated at a point where the ParabolicQfactor was defined correctly.
    fn instantiation_with_psi_last() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        assert_abs_diff_eq!(qfactor.psi_last(), 0.45);
        assert_relative_eq!(qfactor.psip_last(), 0.25921022097041035, epsilon = 1e-12);

        assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);
    }

    #[test]
    /// Values calculated at a point where the ParabolicQfactor was defined correctly.
    fn instantiation_with_psip_last() {
        let qfactor = ParabolicQfactor::new(
            1.1,
            3.9,
            LastClosedFluxSurface::Poloidal(0.25921022097041035),
        );
        assert_abs_diff_eq!(qfactor.psip_last(), 0.25921022097041035);
        assert_relative_eq!(qfactor.psi_last(), 0.45, epsilon = 1e-12);

        assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);
    }

    #[test]
    fn inverse_q() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        let acc = &mut Accelerator::new();

        let psi = 0.1;
        let q_of_psi = qfactor.q_of_psi(psi, acc).unwrap();
        let psi_inverse = qfactor.psi_of_q(q_of_psi, acc).unwrap();
        assert_relative_eq!(psi, psi_inverse, epsilon = 1e-12);

        let psip = 0.1;
        let q_of_psip = qfactor.q_of_psip(psip, acc).unwrap();
        let psip_inverse = qfactor.psip_of_q(q_of_psip, acc).unwrap();
        assert_relative_eq!(psip, psip_inverse, epsilon = 1e-12);
    }
}

#[cfg(test)]
mod test_derivatives_closeness {
    use crate::extract::TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn unity_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Poloidal(0.45));
        test_dpsi_dpsip_q_closeness(&qfactor);
    }

    #[test]
    fn unity_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.45));
        test_dpsip_dpsi_iota_closeness(&qfactor);
    }

    #[test]
    fn parabolic_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        test_dpsi_dpsip_q_closeness(&qfactor);
    }

    #[test]
    fn parabolic_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        test_dpsip_dpsi_iota_closeness(&qfactor);
    }

    #[test]
    fn nc_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        test_dpsi_dpsip_q_closeness(&qfactor);
    }

    #[test]
    fn nc_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        test_dpsip_dpsi_iota_closeness(&qfactor);
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Bad);

        assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Bad);
        assert!(qfactor.q_of_psi_interp.is_some());
        assert!(qfactor.q_of_psip_interp.is_none());
        assert!(qfactor.psip_of_psi_interp.is_some());
        assert!(qfactor.psi_of_psip_interp.is_none());

        assert!(qfactor.psi_array().is_some());
        assert!(qfactor.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        assert!(qfactor.q_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.iota_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.psip_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.dpsip_dpsi(0.01, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EvalError::UndefinedEvaluation as err;
        matches!(qfactor.q_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.iota_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.psi_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.dpsi_dpsip(0.01, &mut acc), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

        assert_eq!(qfactor.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);
        assert!(qfactor.q_of_psi_interp.is_none());
        assert!(qfactor.q_of_psip_interp.is_some());
        assert!(qfactor.psip_of_psi_interp.is_none());
        assert!(qfactor.psi_of_psip_interp.is_some());

        assert!(qfactor.psi_array().is_some());
        assert!(qfactor.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        assert!(qfactor.q_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.iota_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.psi_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.dpsi_dpsip(0.01, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EvalError::UndefinedEvaluation as err;
        matches!(qfactor.q_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.iota_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.psip_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.dpsip_dpsi(0.01, &mut acc), Err(err(..)));
    }
}
