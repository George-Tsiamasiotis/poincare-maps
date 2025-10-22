use crate::Harmonic;
use crate::Result;

use rsl_interpolation::Accelerator;

/// A sum of different perturbation harmonics.
pub struct Perturbation {
    harmonics: Vec<Harmonic>,
}

/// Creation and 'data extraction'
impl Perturbation {
    pub fn from_harmonics(harmonics: Vec<Harmonic>) -> Self {
        Self { harmonics }
    }

    pub fn get_harmonics(&self) -> Vec<Harmonic> {
        self.harmonics.clone()
    }
}

/// Interpolation
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
    /// let path = PathBuf::from("../data.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 1.0, 0.0)?,
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?,
    /// ];
    /// let per = Perturbation::from_harmonics(harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let p = per.p(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn p(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        self.harmonics.iter().try_fold(0.0, |p, harmonic| {
            harmonic.h(psip, theta, zeta, acc).map(|v| p + v)
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
    /// let path = PathBuf::from("../data.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 1.0, 0.0)?,
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?,
    /// ];
    /// let per = Perturbation::from_harmonics(harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let dp_dpsip = per.dp_dpsip(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dpsip(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        self.harmonics.iter().try_fold(0.0, |p, harmonic| {
            harmonic.dh_dpsip(psip, theta, zeta, acc).map(|v| p + v)
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
    /// let path = PathBuf::from("../data.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 1.0, 0.0)?,
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?,
    /// ];
    /// let per = Perturbation::from_harmonics(harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let dp_dtheta = per.dp_dtheta(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        self.harmonics.iter().try_fold(0.0, |p, harmonic| {
            harmonic.dh_dtheta(psip, theta, zeta, acc).map(|v| p + v)
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
    /// let path = PathBuf::from("../data.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 1.0, 0.0)?,
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?,
    /// ];
    /// let per = Perturbation::from_harmonics(harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let dp_dzeta = per.dp_dzeta(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dzeta(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        self.harmonics.iter().try_fold(0.0, |p, harmonic| {
            harmonic.dh_dzeta(psip, theta, zeta, acc).map(|v| p + v)
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
    /// let path = PathBuf::from("../data.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 1.0, 0.0)?,
    ///     Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?,
    /// ];
    /// let per = Perturbation::from_harmonics(harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let dp_dt = per.dp_dt(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dp_dt(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        self.harmonics.iter().try_fold(0.0, |p, harmonic| {
            harmonic.dh_dt(psip, theta, zeta, acc).map(|v| p + v)
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_summation() {
        let path = PathBuf::from("../data.nc");
        let per1 = Perturbation::from_harmonics(vec![
            Harmonic::from_dataset(&path, "akima", 2.0, 1.0, 0.0).unwrap(),
        ]);
        let per2 = Perturbation::from_harmonics(vec![
            Harmonic::from_dataset(&path, "akima", 2.0, 1.0, 0.0).unwrap(),
            Harmonic::from_dataset(&path, "akima", 2.0, 1.0, 0.0).unwrap(),
            Harmonic::from_dataset(&path, "akima", 2.0, 1.0, 0.0).unwrap(),
        ]);

        let mut acc = Accelerator::new();
        let psip = per1.harmonics[0].psip_wall / 2.0;
        let theta = 1.0;
        let zeta = 1.0;

        assert_eq!(
            3.0 * per1.p(psip, theta, zeta, &mut acc).unwrap(),
            per2.p(psip, theta, zeta, &mut acc).unwrap(),
        );
        assert_eq!(
            3.0 * per1.dp_dpsip(psip, theta, zeta, &mut acc).unwrap(),
            per2.dp_dpsip(psip, theta, zeta, &mut acc).unwrap(),
        );
        assert_eq!(
            3.0 * per1.dp_dtheta(psip, theta, zeta, &mut acc).unwrap(),
            per2.dp_dtheta(psip, theta, zeta, &mut acc).unwrap(),
        );
        assert_eq!(
            3.0 * per1.dp_dzeta(psip, theta, zeta, &mut acc).unwrap(),
            per2.dp_dzeta(psip, theta, zeta, &mut acc).unwrap(),
        );
        assert_eq!(
            3.0 * per1.dp_dt(psip, theta, zeta, &mut acc).unwrap(),
            per2.dp_dt(psip, theta, zeta, &mut acc).unwrap(),
        );
    }

    #[test]
    fn test_perturbation_misc() {
        let path = PathBuf::from("../data.nc");
        let harmonics = vec![Harmonic::from_dataset(&path, "akima", 2.0, 1.0, 0.0).unwrap()];
        let per = Perturbation::from_harmonics(harmonics);

        let _ = per.get_harmonics();
    }
}
