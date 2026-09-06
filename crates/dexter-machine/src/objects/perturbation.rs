//! Representation of a total perturbation, a sum of multiple modes.

use crate::{DynMode, DynModeCache, EvalError, MagneticFlux};

/// Reference to a vector of [`DynMode`] objects.
pub type DynModes = Vec<DynMode>;

/// Reference to a vector of [`DynModeCache`] objects.
pub type DynModeCaches = Vec<DynModeCache>;

/// A sum of an arbitrary number of [`Modes`](crate::Mode).
///
/// It has the general form `Σ{ Φ(n,m)(ψ/ψp, θ, ζ, t)}`, where `Φ(n,m)` the different modes.
///
/// Evaluation methods return the sum of the every corresponding evaluation method in each mode.
///
/// The modes are stored in the same order as passed on [`Perturbation::new`].
#[non_exhaustive]
pub struct Perturbation(DynModes);

impl Perturbation {
    /// Returns a Perturbation without any modes, corresponding to an unperturbed state.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let perturbation = Perturbation::zero();
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub const fn zero() -> Self {
        Self(vec![])
    }
}

impl Perturbation {
    /// Creates a new Perturbation from an arbitrary number of modes of any type.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// // from analytical flute modes
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    ///
    /// // from numerical flute modes
    /// let path = PathBuf::from("./netcdf.nc");
    /// let typ = Interpolation1dType::Cubic;
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(NcFluteModeBuilder::new(&path, typ, 2, 1).build()?),
    ///     Box::new(NcFluteModeBuilder::new(&path, typ, 2, 2).build()?),
    ///     Box::new(NcFluteModeBuilder::new(&path, typ, 3, 1).build()?),
    ///     Box::new(NcFluteModeBuilder::new(&path, typ, 3, 2).build()?),
    /// ]);
    ///
    /// // Or any combination of the two
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(NcFluteModeBuilder::new(&path, typ, 2, 1).build()?),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 19, 30, 0.0)),
    /// ]);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn new(modes: Vec<DynMode>) -> Self {
        Self(modes)
    }

    /// Returns the number of the contained modes.
    #[must_use]
    pub fn count(&self) -> usize {
        self.0.len()
    }

    /// Generates a [`Vec`] of the corresponding caching objects.
    ///
    /// The vec has the same length as the `modes` vec, and should be used for evaluations.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let mut caches: DynModeCaches = perturbation.generate_caches();
    /// assert_eq!(perturbation.count(), caches.len());
    ///
    /// let psi = MagneticFlux::Toroidal(0.01);
    /// let p_of_spi = perturbation.eval_p(psi, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn generate_caches(&self) -> DynModeCaches {
        self.0.iter().map(|mode| mode.generate_cache()).collect()
    }
}

/// Evaluations.
impl Perturbation {
    /// Calculates the Perturbation's value as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let psi = MagneticFlux::Toroidal(0.01);
    ///
    /// let p_of_psi = perturbation.eval_p(psi, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the evaluations fail for any reason.
    pub fn eval_p(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut DynModeCaches,
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, mode) = tuple;
                mode.eval_m(flux, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `ψ`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let dp_dpsip = perturbation.eval_deriv_flux(psip, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the evaluations fail for any reason.
    pub fn eval_deriv_flux(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut DynModeCaches,
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, mode) = tuple;
                mode.eval_deriv_flux(flux, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `θ`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let dp_theta = perturbation.eval_deriv_theta(psip, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the evaluations fail for any reason.
    pub fn eval_deriv_theta(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut DynModeCaches,
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, mode) = tuple;
                mode.eval_deriv_theta(flux, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `ζ`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let psi = MagneticFlux::Toroidal(0.01);
    ///
    /// let dp_dzetta = perturbation.eval_deriv_zeta(psi, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the evaluations fail for any reason.
    pub fn eval_deriv_zeta(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut DynModeCaches,
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, mode) = tuple;
                mode.eval_deriv_zeta(flux, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `t`, as a function of `(ψp, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let dp_dt = perturbation.eval_deriv_t(psip, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the evaluations fail for any reason.
    pub fn eval_deriv_t(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut DynModeCaches,
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, modes) = tuple;
                modes
                    .eval_deriv_t(flux, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }
}

impl<Idx> std::ops::Index<Idx> for Perturbation
where
    Idx: std::slice::SliceIndex<[DynMode], Output = DynMode>,
{
    type Output = DynMode;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.0[index]
    }
}

impl std::fmt::Debug for Perturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Perturbation")
            .field("number of modes", &self.0.len())
            .finish()
    }
}

#[cfg(test)]
mod perturbation_evals {
    use std::path::PathBuf;

    use crate::MagneticFlux::Toroidal;
    use crate::extract::TEST_NETCDF_PATH;
    use crate::*;
    use approx::assert_relative_eq;

    use super::*;

    fn create_flute_mode_perturbation() -> Perturbation {
        let lcfs = LastClosedFluxSurface::Toroidal(0.45);
        Perturbation::new(vec![
            Box::new(FluteMode::new(1.0, lcfs, 2, 3, 4.0)),
            Box::new(FluteMode::new(5.0, lcfs, 6, 7, 8.0)),
            Box::new(FluteMode::new(9.0, lcfs, 1, 2, 3.0)),
        ])
    }

    #[rustfmt::skip]
    fn create_nc_flute_mode_perturbation() -> Perturbation {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let typ = Interpolation1dType::Cubic;
        Perturbation::new(vec![
            Box::new(NcFluteModeBuilder::new(&path, typ, 2, 1).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, typ, 2, 1).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, typ, 2, 1).build().unwrap()),
        ])
    }

    #[test]
    fn empty_perturbation() {
        let per = Perturbation::zero();
        let c = &mut per.generate_caches();
        let (p, t) = (Toroidal(0.01), 0.0); // Not used
        assert_eq!(per.eval_p(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.eval_deriv_flux(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.eval_deriv_theta(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.eval_deriv_zeta(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.eval_deriv_t(p, 10.0, 20.0, t, c).unwrap(), 0.0);
    }

    #[test]
    #[rustfmt::skip]
    fn cos_perturbation_evals() {
        let per = create_flute_mode_perturbation();
        let c = &mut per.generate_caches();
        let (p, theta, zeta, t) = (Toroidal(0.1), 0.2, 0.3, 0.0); // Not used
        let eps = 1e-10;
        assert_relative_eq!(per.eval_p(p, theta, zeta, t, c).unwrap(), -2.46342903902, epsilon = eps);
        assert_relative_eq!(per.eval_deriv_flux(p, theta, zeta, t, c).unwrap(), -12.3171451951, epsilon = eps);
        assert_relative_eq!(per.eval_deriv_theta(p, theta, zeta, t, c).unwrap(), -12.1655445266, epsilon = eps);
        assert_relative_eq!(per.eval_deriv_zeta(p, theta, zeta, t, c).unwrap(), 15.9054673268, epsilon = eps);
    }

    #[test]
    #[rustfmt::skip]
    fn nc_perturbation_evals() {
        let per = create_nc_flute_mode_perturbation();
        let c = &mut per.generate_caches();
        let (p, t) = (Toroidal(0.01), 0.0); // Not used
        assert!(per.eval_p(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.eval_deriv_flux(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.eval_deriv_theta(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.eval_deriv_zeta(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.eval_deriv_t(p, 10.0, 20.0, t, c).unwrap().is_finite());
    }

    #[test]
    #[allow(unused_results)]
    fn perturbation_cache() {
        let per = create_flute_mode_perturbation();
        let mut c = per.generate_caches();
        let (p, t) = (Toroidal(0.01), 0.0); // Not used

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 0);
            assert_eq!(cache.misses(), 0);
        });

        per.eval_p(p, 10.0, 20.0, t, &mut c).unwrap();

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 0);
            assert_eq!(cache.misses(), 1);
        });

        per.eval_p(p, 10.0, 20.0, t, &mut c).unwrap();

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 1);
            assert_eq!(cache.misses(), 1);
        });

        per.eval_deriv_theta(p, 10.0, 20.0, t, &mut c).unwrap();
        per.eval_deriv_theta(p, 10.0, 20.0, t, &mut c).unwrap();
        per.eval_deriv_zeta(p, 10.0, 20.0, t, &mut c).unwrap();
        per.eval_deriv_zeta(p, 10.0, 20.0, t, &mut c).unwrap();

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 5);
            assert_eq!(cache.misses(), 1);
        });
    }
}
