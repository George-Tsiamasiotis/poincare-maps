use std::f64::consts::TAU;
use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline, make_spline};

use crate::Result;

use ndarray::Array1;

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
    /// The initial phase of the harmonic.
    pub phase: f64,

    /// The maximum value of the `α` values.
    pub amax: f64,
    /// The value of the poloidal angle ψp at the wall.
    pub psip_wall: f64,
}

/// Holds the Harmonic's values evalutated at a specific spot. Since all the harmonic's methods are
/// called consecutively over the same coordinates, most terms do not need to be calculated every
/// time.
///
/// Similar to the Accelerators, they are stored inside State, and do not affect the behavior of the
/// equilibrium objects themselves.
///
/// The cache should be cloned in each new state calculated from the Solver.
#[derive(Clone, Default)]
pub struct HarmonicCache {
    pub hits: usize,
    pub misses: usize,
    pub psip: f64,
    pub theta: f64,
    pub zeta: f64,
    pub alpha: f64,
    pub dalpha: f64,
    pub sin: f64,
    pub cos: f64,
}

impl HarmonicCache {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_updated(&mut self, psip: f64, theta: f64, zeta: f64) -> bool {
        if (self.psip == psip) & (self.theta == theta) & (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    pub fn update(
        &mut self,
        h: &Harmonic,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
    ) -> Result<()> {
        // Comparing floats is OK here since they are simply copied between every call.
        self.psip = psip;
        self.theta = theta;
        self.zeta = zeta;
        self.alpha = h.a_spline.eval(psip, acc)?;
        self.dalpha = h.a_spline.eval_deriv(psip, acc)?;
        let mod_arg = (h.m * theta - h.n * zeta) % TAU;
        (self.sin, self.cos) = mod_arg.sin_cos();
        Ok(())
    }
}

/// Creation
impl Harmonic {
    /// Constructs a [`Harmonic`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// The spline is only over the amplitude `α`, of the perturbation, and the rest of the
    /// exrpession is analytic.
    ///
    /// # Note
    ///
    /// The value `ψ = 0.0` is prepended at the ψ data array, and the first value of the q array is
    /// prepended (duplicated) in the q array, to assure correct interpolation near the magnetic axis.
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
    pub fn from_dataset(path: &PathBuf, typ: &str, m: f64, n: f64, phase: f64) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        // Make path absolute. Just unwrap, Equilibrium checks if it exists.
        let path = std::path::absolute(path).unwrap();

        let eq = Equilibrium::from_file(&path)?;

        let psip_data = extract_1d_var(&eq.file, PSIP_COORD)?
            .as_standard_layout()
            .to_vec();
        let psip_wall = psip_data.last().copied().unwrap();

        // TODO: update
        let mut a_data = Array1::zeros(psip_data.len());
        for i in 0..psip_data.len() {
            a_data[i] = gaussian(psip_data[i], psip_wall)
        }
        let amax = *a_data
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&0.0);

        let a_spline = make_spline(typ, &psip_data, &a_data.to_vec())?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            a_spline,
            psip_wall,
            m,
            n,
            phase: phase % TAU,
            amax,
        })
    }
}

/// Interpolation
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
        psip: f64,
        theta: f64,
        zeta: f64,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(&self, psip, theta, zeta, acc)?
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
        psip: f64,
        theta: f64,
        zeta: f64,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(&self, psip, theta, zeta, acc)?
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
        psip: f64,
        theta: f64,
        zeta: f64,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(&self, psip, theta, zeta, acc)?
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
        psip: f64,
        theta: f64,
        zeta: f64,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(&self, psip, theta, zeta, acc)?
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
        psip: f64,
        theta: f64,
        zeta: f64,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        // Time-independent perturbations at the moment.
        Ok(0.0)
    }
}

/// Data extraction
impl Harmonic {
    /// Returns the `psip` coordinate data as a 1D array.
    pub fn psip_data(&self) -> Array1<f64> {
        Array1::from_vec(self.a_spline.xa.to_vec())
    }

    /// Returns the `α` perturbation amplitude data as a 1D array.
    pub fn a_data(&self) -> Array1<f64> {
        Array1::from_vec(self.a_spline.ya.to_vec())
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
            amax: self.amax,
            psip_wall: self.psip_wall,
        }
    }
}

impl std::fmt::Debug for Harmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Harmonic")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall))
            .field("m", &self.m)
            .field("n", &self.n)
            .field("φ", &self.phase)
            .field("α_max", &self.amax)
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

        assert_eq!(harmonic.psip_data().ndim(), 1);
        assert_eq!(harmonic.a_data().ndim(), 1);
    }

    #[test]
    fn test_harmonic_misc() {
        let path = PathBuf::from("../data.nc");
        let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0).unwrap();
        let _ = harmonic.clone();
        let _ = format!("{harmonic:?}");
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
