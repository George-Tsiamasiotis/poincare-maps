use crate::Result;
use crate::equilibrium::Harmonic;

use rsl_interpolation::Accelerator;

use pyo3::prelude::*;
use pyo3::types::PyList;

/// A sum of different perturbation harmonics.
#[pyclass(frozen, immutable_type)]
pub struct Perturbation {
    harmonics: Vec<Harmonic>,
}

#[pymethods]
impl Perturbation {
    /// Creates a new [`Perturbation`]
    ///
    /// Wrapper around [`Perturbation::from_dataset`]. This is a workaround to return a [`PyErr`].
    #[coverage(off)]
    #[new]
    pub fn new_py<'py>(harmonics: Bound<'py, PyList>) -> Self {
        let harmonics_vec: Vec<Harmonic> = harmonics.iter().map(|h| h.extract().unwrap()).collect();
        dbg!(harmonics_vec.len());
        Self::from_harmonics(harmonics_vec)
    }
}

impl Perturbation {
    pub fn from_harmonics(harmonics: Vec<Harmonic>) -> Self {
        Self {
            harmonics: harmonics,
        }
    }
}

/// Evaluation functions
impl Perturbation {
    /// Calculates the Perturbation `Σ{ α(n,m)(ψp) * cos(mθ-nζ+φ0) }`.
    ///
    /// TODO: EXAMPLE
    pub fn p(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let mut p = 0.0;
        self.harmonics.iter().try_fold(0.0, |_, harmonic| {
            harmonic.h(psip, theta, zeta, acc).inspect(|val| p += val)
        })?;
        Ok(p)
    }

    /// Calculates the Perturbation's derivative with respect to `ψp`.
    ///
    /// TODO: EXAMPLE
    pub fn dp_dpsip(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let mut p = 0.0;
        self.harmonics.iter().try_fold(0.0, |_, harmonic| {
            harmonic
                .dh_dpsip(psip, theta, zeta, acc)
                .inspect(|val| p += val)
        })?;
        Ok(p)
    }

    /// Calculates the Perturbation's derivative with respect to `θ`.
    ///
    /// TODO: EXAMPLE
    pub fn dp_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        let mut p = 0.0;
        self.harmonics.iter().try_fold(0.0, |_, harmonic| {
            harmonic
                .dh_dtheta(psip, theta, zeta, acc)
                .inspect(|val| p += val)
        })?;
        Ok(p)
    }

    /// Calculates the Perturbation's derivative with respect to `ζ`.
    ///
    /// TODO: EXAMPLE
    pub fn dp_dzeta(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let mut p = 0.0;
        self.harmonics.iter().try_fold(0.0, |_, harmonic| {
            harmonic
                .dh_dzeta(psip, theta, zeta, acc)
                .inspect(|val| p += val)
        })?;
        Ok(p)
    }

    /// Calculates the Perturbation's derivative with respect to `t`.
    ///
    /// TODO: EXAMPLE
    pub fn dp_dt(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let mut p = 0.0;
        self.harmonics.iter().try_fold(0.0, |_, harmonic| {
            harmonic
                .dh_dt(psip, theta, zeta, acc)
                .inspect(|val| p += val)
        })?;
        Ok(p)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_summation() {
        let path = PathBuf::from("./data.nc");
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

        dbg!(per1.p(psip, theta, zeta, &mut acc).unwrap(),);
        dbg!(per2.p(psip, theta, zeta, &mut acc).unwrap(),);
        // assert_eq!(
        //     2.0 * per1.p(psip, theta, zeta, &mut acc).unwrap(),
        //     per2.p(psip, theta, zeta, &mut acc).unwrap(),
        // );
    }
}
