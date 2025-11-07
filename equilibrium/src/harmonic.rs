use std::f64::consts::TAU;
use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline, make_spline};
use utils::array1D_getter_impl;

use crate::Result;
use crate::{Flux, Radians};

use ndarray::Array1;
use safe_unwrap::safe_unwrap;

/// Single perturbation harmonic reconstructed from a netCDF file.
///
/// The harmonic has the form of `α(ψp) * cos(mθ-nζ+φ0)`, where `α(ψp)` is calculated by
/// interpolation over some numerical data.
pub struct Harmonic {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// Interpolation type.
    pub typ: String,
    /// Spline over the perturbation amplitude `α` data, as a function of ψp.
    pub a_spline: DynSpline<f64>,
    /// The `θ` frequency number.
    pub m: f64,
    /// The `ζ` frequency number.
    pub n: f64,
    /// The phase offset of the harmonic.
    pub phase: Radians,
}

/// Holds the Harmonic's values evalutated at a specific point.
///
/// Since all the harmonic's methods are called consecutively over the same coordinates, most terms
/// do not need to be calculated every time.
///
/// Similar to the Accelerators, they are stored inside State, and do not affect the behavior of the
/// equilibrium objects themselves.
///
/// The cache should be cloned in each new state calculated from the Solver.
#[derive(Clone, Default)]
pub struct HarmonicCache {
    pub hits: usize,
    pub misses: usize,
    pub psip: Flux,
    pub theta: Radians,
    pub zeta: Radians,
    pub alpha: f64,
    pub dalpha: f64,
    pub sin: f64,
    pub cos: f64,
}

impl HarmonicCache {
    pub fn new() -> Self {
        Self::default()
    }

    /// Checks if the cache's fields are valid.
    ///
    /// Comparing floats is OK here since they are simply copied between every call, and we want
    /// the check to fail with the slightest difference.
    pub fn is_updated(&mut self, psip: Flux, theta: Radians, zeta: Radians) -> bool {
        if (self.psip == psip) & (self.theta == theta) & (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    /// Updates the cache's fields.
    pub fn update(
        &mut self,
        h: &Harmonic,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        acc: &mut Accelerator,
    ) -> Result<()> {
        self.psip = psip;
        self.theta = theta;
        self.zeta = zeta;
        self.alpha = h.a_spline.eval(psip, acc)?;
        self.dalpha = h.a_spline.eval_deriv(psip, acc)?;
        let mod_arg = (h.m * theta - h.n * zeta + h.phase) % TAU;
        (self.sin, self.cos) = mod_arg.sin_cos();
        Ok(())
    }
}

// Creation
impl Harmonic {
    /// Constructs a [`Harmonic`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// The spline is only over the amplitude `α`, of the perturbation, and the rest of the
    /// exrpession is analytic.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str, m: f64, n: f64, phase: Radians) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;

        let eq = Equilibrium::from_file(&path)?;

        let psip_data = extract_1d_var(&eq.file, PSIP_COORD)?
            .as_standard_layout()
            .to_owned();

        // TODO: update
        let mut a_data = Array1::zeros(psip_data.len());
        for i in 0..psip_data.len() {
            a_data[i] = gaussian(psip_data[i], psip_data.last().copied().unwrap())
        }

        let a_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", a_data.as_slice()),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            a_spline,
            m,
            n,
            phase: phase % TAU,
        })
    }
}

// Interpolation
impl Harmonic {
    /// Calculates the harmonic `a(ψp) * <analytical term>`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let h = harmonic.h(0.015, 2.0*PI, 0.0, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn h(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * cache.cos)
    }

    /// Calculates the harmonic derivative `𝜕h/𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dpsip = harmonic.dh_dpsip(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dh_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.dalpha * cache.cos)
    }

    /// Calculates the harmonic derivative `𝜕h/𝜕θ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dtheta = harmonic.dh_dtheta(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dh_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * (-self.m) * cache.sin)
    }

    /// Calculates the perturbation derivative `𝜕h/𝜕ζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dzeta = harmonic.dh_dzeta(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dh_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * self.n * cache.sin)
    }

    /// Calculates the perturbation derivative `𝜕h/𝜕t`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dt = harmonic.dh_dt(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    #[allow(unused_variables)]
    pub fn dh_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        // Time-independent perturbations at the moment.
        Ok(0.0)
    }
}

// Data extraction
impl Harmonic {
    array1D_getter_impl!(psip_data, a_spline.xa, Flux);
    array1D_getter_impl!(a_data, a_spline.ya, Flux);

    /// Returns the value of the poloidal angle ψp at the wall.
    pub fn psip_wall(&self) -> Flux {
        safe_unwrap!("ya is non-empty", self.a_spline.xa.last().copied())
    }

    /// Returns the value of the `m` mode number.
    pub fn m(&self) -> f64 {
        self.m
    }

    /// Returns the value of the `n` mode number.
    pub fn n(&self) -> f64 {
        self.n
    }

    /// Returns the value of the harmonic's phase offset.
    pub fn phase(&self) -> f64 {
        self.phase
    }
}

/// A simple gaussian distribution to emulate reconstructed perturbations.
/// TODO: remove
fn gaussian(psip: f64, psip_wall: f64) -> f64 {
    use std::f64::consts::TAU;

    let scale = 2e-2;
    let mu = psip_wall / 2.0;
    let sigma = psip_wall / 4.0;

    scale * (1.0 / (TAU * sigma).sqrt()) * (-(psip - mu).powi(2) / (2.0 * sigma.powi(2))).exp()
}

impl Clone for Harmonic {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            typ: self.typ.clone(),
            a_spline: make_spline(&self.typ, &self.a_spline.xa, &self.a_spline.ya)
                .expect("Could not clone spline."),
            m: self.m,
            n: self.n,
            phase: self.phase,
        }
    }
}

impl std::fmt::Debug for Harmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Harmonic")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall()))
            .field("m", &self.m)
            .field("n", &self.n)
            .field("φ", &self.phase)
            .finish()
    }
}
impl std::fmt::Debug for HarmonicCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HarmonicCache")
            .field("hits  ", &self.hits)
            .field("misses", &self.misses)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_data_extraction() {
        let path = PathBuf::from("../data.nc");
        let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0).unwrap();
        let _: f64 = harmonic.m();
        let _: f64 = harmonic.n();
        let _: f64 = harmonic.phase();

        assert_eq!(harmonic.psip_data().ndim(), 1);
        assert_eq!(harmonic.a_data().ndim(), 1);
    }

    #[test]
    fn test_cache_update() {
        let path = PathBuf::from("../data.nc");
        let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0).unwrap();
        let mut acc = Accelerator::new();
        let mut cache = HarmonicCache::new();

        let (psip, theta, zeta) = (0.015, 0.0, 3.14);
        harmonic.h(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        harmonic
            .dh_dpsip(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        harmonic
            .dh_dtheta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        harmonic
            .dh_dzeta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        harmonic
            .dh_dt(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        let (psip, theta, zeta) = (0.01, 0.01, 3.15);
        harmonic.h(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        harmonic
            .dh_dpsip(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        harmonic
            .dh_dtheta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        harmonic
            .dh_dzeta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        harmonic
            .dh_dt(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
    }

    #[test]
    fn test_harmonic_misc() {
        let path = PathBuf::from("../data.nc");
        let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0).unwrap();
        let _ = harmonic.clone();
        let _ = format!("{harmonic:?}");
    }
}
