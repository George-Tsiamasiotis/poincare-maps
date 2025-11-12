use crate::Harmonic;
use crate::HarmonicCache;

use crate::Result;
use crate::{Flux, Radians};

use rsl_interpolation::Accelerator;

/// A sum of different perturbation harmonics.
pub struct Perturbation {
    pub harmonics: Vec<Harmonic>,
}

// Creation and data extraction
impl Perturbation {
    pub fn from_harmonics(harmonics: &[Harmonic]) -> Self {
        Self {
            harmonics: harmonics.into(),
        }
    }

    pub fn get_harmonics(&self) -> Vec<Harmonic> {
        self.harmonics.clone()
    }
}

// Interpolation
impl Perturbation {
    /// Calculates the Perturbation `Σ{ α(n,m)(ψp) * cos(mθ-nζ+φ0) }`.
    ///
    /// If an empty Vec is passed, it is equivalent to the non-perturbed system.
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
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let p = per.p(0.015, 2.0*PI, PI, &mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn p(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [HarmonicCache],
        acc: &mut Accelerator,
    ) -> Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .h(psip, theta, zeta, &mut caches[index], acc)
                .map(|v| p + v)
        })
    }

    /// Calculates the Perturbation's derivative with respect to `ψp`,
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
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dpsip = per.dp_dpsip(0.015, 2.0*PI, PI,&mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut [HarmonicCache],
        acc: &mut Accelerator,
    ) -> Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dpsip(psip, theta, zeta, &mut cache[index], acc)
                .map(|v| p + v)
        })
    }

    /// Calculates the Perturbation's derivative with respect to `θ`.
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
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dtheta = per.dp_dtheta(0.015, 2.0*PI, PI, &mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut [HarmonicCache],
        acc: &mut Accelerator,
    ) -> Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dtheta(psip, theta, zeta, &mut cache[index], acc)
                .map(|v| p + v)
        })
    }

    /// Calculates the Perturbation's derivative with respect to `ζ`.
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
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dzeta = per.dp_dzeta(0.015, 2.0*PI, PI, &mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut [HarmonicCache],
        acc: &mut Accelerator,
    ) -> Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dzeta(psip, theta, zeta, &mut cache[index], acc)
                .map(|v| p + v)
        })
    }

    /// Calculates the Perturbation's derivative with respect to `t`.
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
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dt = per.dp_dt(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut [HarmonicCache],
        acc: &mut Accelerator,
    ) -> Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dt(psip, theta, zeta, &mut cache[index], acc)
                .map(|v| p + v)
        })
    }
}

impl std::fmt::Debug for Perturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for harmonic in self.get_harmonics() {
            let _ = harmonic.fmt(f);
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use config::STUB_NETCDF_PATH;

    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_summation() {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        let per1 = Perturbation::from_harmonics(&vec![
            Harmonic::from_dataset(&path, "akima", 2, 1).unwrap(),
        ]);
        let per2 = Perturbation::from_harmonics(&vec![
            Harmonic::from_dataset(&path, "akima", 2, 1).unwrap(),
            Harmonic::from_dataset(&path, "akima", 2, 1).unwrap(),
            Harmonic::from_dataset(&path, "akima", 2, 1).unwrap(),
        ]);

        let mut acc = Accelerator::new();
        // Normally, this would happen inside State.
        let mut hcache1 = vec![HarmonicCache::default(); per1.harmonics.len()];
        let mut hcache2 = vec![HarmonicCache::default(); per2.harmonics.len()];
        let psip = per1.harmonics[0].psip_wall() / 2.0;
        let theta = 1.0;
        let zeta = 1.0;

        assert_eq!(
            3.0 * per1.p(psip, theta, zeta, &mut hcache1, &mut acc).unwrap(),
            per2.p(psip, theta, zeta, &mut hcache2, &mut acc).unwrap(),
        );
        assert_eq!(
            3.0 * per1
                .dp_dpsip(psip, theta, zeta, &mut hcache1, &mut acc)
                .unwrap(),
            per2.dp_dpsip(psip, theta, zeta, &mut hcache2, &mut acc)
                .unwrap(),
        );
        assert_eq!(
            3.0 * per1
                .dp_dtheta(psip, theta, zeta, &mut hcache1, &mut acc)
                .unwrap(),
            per2.dp_dtheta(psip, theta, zeta, &mut hcache2, &mut acc)
                .unwrap(),
        );
        assert_eq!(
            3.0 * per1
                .dp_dzeta(psip, theta, zeta, &mut hcache1, &mut acc)
                .unwrap(),
            per2.dp_dzeta(psip, theta, zeta, &mut hcache2, &mut acc)
                .unwrap(),
        );
        assert_eq!(
            3.0 * per1
                .dp_dt(psip, theta, zeta, &mut hcache1, &mut acc)
                .unwrap(),
            per2.dp_dt(psip, theta, zeta, &mut hcache2, &mut acc)
                .unwrap(),
        );
    }

    #[test]
    fn test_perturbation_misc() {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        let harmonics = vec![Harmonic::from_dataset(&path, "akima", 2, 1).unwrap()];
        let per = Perturbation::from_harmonics(&harmonics);

        let _ = per.get_harmonics();
    }
}
