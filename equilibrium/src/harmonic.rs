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
/// The harmonic has the form of `α(ψp) * cos(mθ-nζ+φ(ψp))`, where `α(ψp)` and `φ(ψp)` are calculated by
/// interpolation over numerical data.
pub struct Harmonic {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    pub typ: String,
    /// Spline over the perturbation amplitude `α` data, as a function of ψp.
    pub a_spline: DynSpline<f64>,
    /// Spline over the perturbation amplitude `φ` data, as a function of ψp.
    pub phase_spline: DynSpline<f64>,
    /// The `θ` frequency number.
    pub m: i64,
    /// The `ζ` frequency number.
    pub n: i64,

    // Used in the actual calculations
    _m: f64,
    _n: f64,
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
    pub phase: f64,
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
        self.phase = h.phase_spline.eval(psip, acc)?;
        self.dalpha = h.a_spline.eval_deriv(psip, acc)?;
        let mod_arg = (h._m * self.theta - h._n * self.zeta + self.phase) % TAU;
        (self.sin, self.cos) = mod_arg.sin_cos();
        Ok(())
    }
}

// Creation
impl Harmonic {
    /// Constructs a [`Harmonic`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// The spline is only over the amplitude `α` and phase `φ`, of the perturbation. The rest of the
    /// exrpession is analytic.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str, m: i64, n: i64) -> Result<Self> {
        use crate::extract::*;
        use config::netcdf_fields::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP)?.as_standard_layout().to_owned();
        let (a_data, phase_data) = extract_harmonic_arrays(&f, m, n)?;

        let a_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", a_data.as_slice()),
        )?;
        let phase_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", phase_data.as_slice()),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            a_spline,
            phase_spline,
            m,
            n,
            _m: m as f64,
            _n: n as f64,
        })
    }
}

// Interpolation
impl Harmonic {
    /// Calculates the harmonic `α(ψp) * cos(mθ-nζ+φ(ψp))`.
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
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
        Ok(cache.alpha * (-self._m) * cache.sin)
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
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
        Ok(cache.alpha * self._n * cache.sin)
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
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
    pub fn m(&self) -> i64 {
        self.m
    }

    /// Returns the value of the `n` mode number.
    pub fn n(&self) -> i64 {
        self.n
    }
}

impl Clone for Harmonic {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            typ: self.typ.clone(),
            a_spline: make_spline(&self.typ, &self.a_spline.xa, &self.a_spline.ya)
                .expect("Could not clone spline."),
            phase_spline: make_spline(&self.typ, &self.phase_spline.xa, &self.phase_spline.ya)
                .expect("Could not clone spline."),
            m: self.m,
            n: self.n,
            _m: self._m,
            _n: self._n,
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
    use config::STUB_NETCDF_PATH;

    fn get_test_dataset_path() -> PathBuf {
        PathBuf::from(STUB_NETCDF_PATH)
    }

    #[test]
    fn test_data_extraction() {
        let path = get_test_dataset_path();
        let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2).unwrap();
        let _: i64 = harmonic.m();
        let _: i64 = harmonic.n();

        assert_eq!(harmonic.psip_data().ndim(), 1);
        assert_eq!(harmonic.a_data().ndim(), 1);
    }

    #[test]
    fn test_cache_update() {
        let path = get_test_dataset_path();
        let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2).unwrap();
        let mut acc = Accelerator::new();
        let mut cache = HarmonicCache::new();

        // dh_dt does not update the cache
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
        assert_eq!(cache.misses, 1);
        assert_eq!(cache.hits, 3);
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
        assert_eq!(cache.misses, 2);
        assert_eq!(cache.hits, 6);
    }

    #[test]
    fn test_harmonic_misc() {
        let path = get_test_dataset_path();
        let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2).unwrap();
        let _ = harmonic.clone();
        let _ = format!("{harmonic:?}");
    }
}
