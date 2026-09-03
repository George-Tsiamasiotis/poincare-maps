//! Representation of analytical Flute Modes.

use rsl_interpolation::Accelerator;
use std::f64::consts::TAU;

use crate::{
    DynModeCache, EvalError, FluxCoordinateState, LastClosedFluxSurface, MachineObject,
    MachineType, Mode, ModeCache,
};
use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    mode_cache_getters_impl,
};

// ===============================================================================================

/// A simple analytical flute mode of the form `ε*sqrt(ψ/ψlast)*cos(mθ-nζ+φ)` or
/// `ε*sqrt(ψp/ψplast)*cos(mθ-nζ+φ)`, where `ε` and `φ` are constants and `ψlast/ψplast` the last
/// closed flux surface.
///
/// The square root is necessary since the modes must behave like the square root of the
/// magnetic flux close to the axis.
///
/// Used in pair with [`FluteModeCache`].
#[non_exhaustive]
#[derive(Clone)]
pub struct FluteMode {
    /// The modes's "amplitude" `ε`. Corresponds the value of the amplitude at the last closed
    /// flux surface.
    epsilon: f64,
    /// The modes's poloidal mode number `m`.
    m: i64,
    /// The modes's toroidal mode number `n`.
    n: i64,
    /// The modes's phase.
    phase: f64,
    /// The last closed flux surface.
    lcfs: LastClosedFluxSurface,
    /// The value of the last closed toroidal flux surface, if the mode was defined through the
    /// toroidal flux.
    psi_last: Option<f64>,
    /// The value of the last closed poloidal flux surface, if the mode was defined through the
    /// poloidal flux.
    psip_last: Option<f64>,
}

impl FluteMode {
    /// Creates a new `FluteMode`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// ```
    #[must_use]
    pub fn new(epsilon: f64, lcfs: LastClosedFluxSurface, m: i64, n: i64, phase: f64) -> Self {
        let psi_last: Option<f64>;
        let psip_last: Option<f64>;
        match lcfs {
            LastClosedFluxSurface::Toroidal(last) => {
                psi_last = Some(last);
                psip_last = None;
            }
            LastClosedFluxSurface::Poloidal(last) => {
                psi_last = None;
                psip_last = Some(last)
            }
        }

        Self {
            epsilon,
            m,
            n,
            phase,
            lcfs,
            psi_last,
            psip_last,
        }
    }

    /// Returns the mode's last closed flux surface.
    #[must_use]
    pub fn lcfs(&self) -> LastClosedFluxSurface {
        self.lcfs
    }

    /// Returns the mode's constant "amplitude" `ε`.
    #[must_use]
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    /// Returns the mode's constant phase `φ`.
    #[must_use]
    pub fn phase(&self) -> f64 {
        self.phase
    }
}

impl MachineObject for FluteMode {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        match self.psi_last {
            Some(_) => FluxCoordinateState::Good,
            None => FluxCoordinateState::Bad,
        }
    }

    fn psip_state(&self) -> FluxCoordinateState {
        match self.psip_last {
            Some(_) => FluxCoordinateState::Good,
            None => FluxCoordinateState::Bad,
        }
    }
}

impl std::fmt::Debug for FluteMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FluteMode")
            .field("epsilon", &self.epsilon)
            .field("LCFS", &self.lcfs())
            .field("poloidal number `m`", &self.m)
            .field("toroidal number `n`", &self.n)
            .field("phase", &self.phase)
            .finish()
    }
}

/// Stores a [`FluteMode`]'s constant parameters and cached quantities.
#[derive(Debug, Clone)]
pub struct FluteModeCache {
    /// The number of cache hits.
    hits: usize,
    /// The number of cache misses.
    misses: usize,
    /// Ordered array with the mode's static parameters.
    ///
    /// cache = [
    ///     0 = epsilon,
    ///     1 = m,
    ///     2 = n,
    ///     3 = phase,
    ///     4 = sqrt(LCFS)
    /// ].
    params: [f64; 5],
    /// Ordered array with the modes intermediate cached values.
    ///
    /// cache = [
    ///     0 = flux,
    ///     1 = theta,
    ///     2 = zeta,
    ///     3 = sqrt(flux)
    ///     4 = modarg
    ///     5 = sin,
    ///     6 = cos
    /// ].
    cache: [f64; 7],
}

impl Default for FluteModeCache {
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            params: [f64::NAN; 5],
            cache: [f64::NAN; 7],
        }
    }
}

impl ModeCache for FluteModeCache {
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
        self.cache[3] = flux.sqrt();
        self.cache[4] =
            (self.params[1] * theta - self.params[2] * zeta + self.params[3]).rem_euclid(TAU);
        (self.cache[5], self.cache[6]) = self.cache[4].sin_cos();
    }

    fn acc(&mut self) -> Option<&mut Accelerator> {
        None
    }

    mode_cache_getters_impl!(FluteModeCache);
}

impl Mode for FluteMode {
    fn psi_last(&self) -> Option<f64> {
        self.psi_last
    }

    fn psip_last(&self) -> Option<f64> {
        self.psip_last
    }

    fn m(&self) -> i64 {
        self.m
    }

    fn n(&self) -> i64 {
        self.n
    }

    fn generate_cache(&self) -> DynModeCache {
        let lcfs_root = match self.lcfs {
            LastClosedFluxSurface::Toroidal(last) | LastClosedFluxSurface::Poloidal(last) => {
                last.sqrt()
            }
        };
        Box::new(FluteModeCache {
            params: [
                self.epsilon,
                self.m as f64,
                self.n as f64,
                self.phase,
                lcfs_root,
            ],
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("α(ψ)".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                self.epsilon * cache.cache()[3] / cache.params()[4]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("α(ψp)".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                self.epsilon * cache.cache()[3] / cache.params()[4]
            ))
        }
    }

    fn phase_of_psi(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(self.phase))
    }

    fn phase_of_psip(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut DynModeCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(debug_assert_is_finite!(self.phase))
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("h(ψ)".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.params()[0] * cache.cache()[3] / cache.params()[4] * cache.cache()[6]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("h(ψp)".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.params()[0] * cache.cache()[3] / cache.params()[4] * cache.cache()[6]
            ))
        }
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψ)/dψ".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.params()[0] / (2.0 * cache.params()[4] * cache.cache()[3]) * cache.cache()[6]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψp)/dψp".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.params()[0] / (2.0 * cache.params()[4] * cache.cache()[3]) * cache.cache()[6]
            ))
        }
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψ)/dθ".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                -cache.params()[1] * cache.params()[0] * cache.cache()[3] / cache.params()[4]
                    * cache.cache()[5]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψp)/dθ".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                -cache.params()[1] * cache.params()[0] * cache.cache()[3] / cache.params()[4]
                    * cache.cache()[5]
            ))
        }
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψ)/dζ".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.params()[2] * cache.params()[0] * cache.cache()[3] / cache.params()[4]
                    * cache.cache()[5]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψp)/dζ".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.params()[2] * cache.params()[0] * cache.cache()[3] / cache.params()[4]
                    * cache.cache()[5]
            ))
        }
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
        Ok(debug_assert_is_finite!(0.0_f64))
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
        Ok(debug_assert_is_finite!(0.0_f64))
    }
}

#[cfg(test)]
mod flute_mode_values {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    #[rustfmt::skip]
    fn desmos_values_toroidal_lcfs() -> Result<(), EvalError> {
        let lcfs = LastClosedFluxSurface::Toroidal(0.45);
        let har = dbg!(FluteMode::new(10.0, lcfs, 3, 2, 1.0));
        assert_eq!(har.psi_state(), FluxCoordinateState::Good);
        assert_eq!(har.psip_state(), FluxCoordinateState::Bad);

        let c = &mut har.generate_cache();

        let p = 0.1; // not used
        let theta = 0.2;
        let zeta = 0.3;
        let t = 0.0; // not used

        let eps = 1e-10;
        assert_relative_eq!(har.ampl_of_psi(p, theta, zeta, t, c)?, 4.71404520791, epsilon = eps);
        assert_relative_eq!(har.phase_of_psi(p, theta, zeta, t, c)?, 1.0, epsilon = eps);
        assert_relative_eq!(har.m_of_psi(p, theta, zeta, t, c)?, 2.5470094958, epsilon = eps);
        assert_relative_eq!(har.dm_dpsi(p, theta, zeta, t, c)?, 12.735047479, epsilon = eps);
        assert_relative_eq!(har.dm_of_psi_dtheta(p, theta, zeta, t, c)?, -11.9001967906, epsilon = eps);
        assert_relative_eq!(har.dm_of_psi_dzeta(p, theta, zeta, t, c)?, 7.93346452706, epsilon = eps);

        assert!(har.ampl_of_psip(p, theta, zeta, t,c).is_err());
        assert!(har.m_of_psip(p, theta, zeta, t,c).is_err());
        assert!(har.dm_dpsip(p, theta, zeta, t,c).is_err());
        assert!(har.dm_of_psip_dtheta(p, theta, zeta, t,c).is_err());
        assert!(har.dm_of_psip_dzeta(p, theta, zeta, t,c).is_err());

        assert_eq!(c.misses(), 1);
        assert_eq!(c.hits(), 4);

        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn desmos_values_poloidal_lcfs() -> Result<(), EvalError> {
        let lcfs = LastClosedFluxSurface::Poloidal(0.45);
        let har = dbg!(FluteMode::new(10.0, lcfs, 3, 2, 1.0));
        assert_eq!(har.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(har.psip_state(), FluxCoordinateState::Good);

        let c = &mut har.generate_cache();

        let p = 0.1; // not used
        let theta = 0.2;
        let zeta = 0.3;
        let t = 0.0; // not used

        let eps = 1e-10;
        assert_relative_eq!(har.ampl_of_psip(p, theta, zeta, t, c)?, 4.71404520791, epsilon = eps);
        assert_relative_eq!(har.phase_of_psip(p, theta, zeta, t, c)?, 1.0, epsilon = eps);
        assert_relative_eq!(har.m_of_psip(p, theta, zeta, t, c)?, 2.5470094958, epsilon = eps);
        assert_relative_eq!(har.dm_dpsip(p, theta, zeta, t, c)?, 12.735047479, epsilon = eps);
        assert_relative_eq!(har.dm_of_psip_dtheta(p, theta, zeta, t, c)?, -11.9001967906, epsilon = eps);
        assert_relative_eq!(har.dm_of_psip_dzeta(p, theta, zeta, t, c)?, 7.93346452706, epsilon = eps);

        assert!(har.ampl_of_psi(p, theta, zeta, t,c).is_err());
        assert!(har.m_of_psi(p, theta, zeta, t,c).is_err());
        assert!(har.dm_dpsi(p, theta, zeta, t,c).is_err());
        assert!(har.dm_of_psi_dtheta(p, theta, zeta, t,c).is_err());
        assert!(har.dm_of_psi_dzeta(p, theta, zeta, t,c).is_err());

        assert_eq!(c.misses(), 1);
        assert_eq!(c.hits(), 4);

        Ok(())
    }
}

#[cfg(test)]
#[expect(unused_results)]
mod flute_mode_cache {

    use super::*;

    #[test]
    fn counts() {
        let lcfs = LastClosedFluxSurface::Toroidal(0.45);
        let mode = dbg!(FluteMode::new(10.0, lcfs, 3, 2, 1.0));
        let c = &mut mode.generate_cache();
        let psi = 0.01; // not checked
        let t = 0.0; // not checked

        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 0);

        mode.phase_of_psi(psi, 0.1, 0.1, t, c).unwrap(); // Does not check
        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 0);

        mode.dm_of_psi_dtheta(0.01, 0.1, 0.1, t, c).unwrap(); // First check
        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 1);

        let theta = 3.14;
        let zeta = 1.0;
        mode.m_of_psi(psi, theta, zeta, t, c).unwrap();
        mode.dm_of_psi_dtheta(psi, theta, zeta, t, c).unwrap();
        mode.dm_of_psi_dzeta(psi, theta, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 2);
        assert_eq!(c.misses(), 2);

        mode.dm_of_psi_dzeta(psi, theta / 2.0, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 2);
        assert_eq!(c.misses(), 3);
    }
}
