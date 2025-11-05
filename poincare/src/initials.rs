use ndarray::Array1;

use particle::{InitialConditions, Particle};
use utils::array1D_getter_impl;

use crate::{PoincareError, Result};

/// Stores the initial conditions arrays.
#[derive(Clone)]
#[non_exhaustive]
pub struct PoincareInit {
    pub thetas: Array1<f64>,
    pub psips: Array1<f64>,
    pub rhos: Array1<f64>,
    pub zetas: Array1<f64>,
    pub mus: Array1<f64>,
}

impl PoincareInit {
    /// Creates a new [`PoincareInit`].
    ///
    /// # Error
    ///
    /// Returns [`PoincareError`] if the input arrays are not of the same length.
    ///
    /// # Example
    /// ```
    /// # use poincare::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let p = PoincareInit::build(
    ///     &[0.0, 0.1, 0.2],
    ///     &[0.0, 0.15, 0.3],
    ///     &[1e-3, 2e-3, 3e-3],
    ///     &[0.0, 0.1, 0.2],
    ///     &[0.0, 0.0, 0.0],
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn build(
        thetas: &[f64],
        psips: &[f64],
        rhos: &[f64],
        zetas: &[f64],
        mus: &[f64],
    ) -> Result<Self> {
        let len = thetas.len();
        if !(thetas.len() == len
            && psips.len() == len
            && rhos.len() == len
            && zetas.len() == len
            && mus.len() == len)
        {
            return Err(PoincareError::InitMismatch);
        }

        Ok(Self {
            thetas: Array1::from_vec(thetas.to_vec()),
            psips: Array1::from_vec(psips.to_vec()),
            rhos: Array1::from_vec(rhos.to_vec()),
            zetas: Array1::from_vec(zetas.to_vec()),
            mus: Array1::from_vec(mus.to_vec()),
        })
    }

    /// Creates a vector with [`Particle`]s initialized over all the initial conditions sets.
    ///
    /// # Example
    /// ```
    /// # use poincare::*;
    /// # use particle::Particle;
    /// #
    /// # fn main() -> Result<()> {
    /// let p = PoincareInit::build(
    ///     &[0.0, 0.1, 0.2],
    ///     &[0.0, 0.15, 0.3],
    ///     &[1e-3, 2e-3, 3e-3],
    ///     &[0.0, 0.1, 0.2],
    ///     &[0.0, 0.0, 0.0],
    /// )?;
    /// let particles: Vec<Particle> = p.to_particles();
    /// # Ok(())
    /// # }
    /// ```
    pub fn to_particles(&self) -> Vec<Particle> {
        let mut particles = Vec::with_capacity(self.len());
        for index in 0..self.len() {
            particles.push(self.particle_from_index(index));
        }
        particles
    }

    /// Creates an [`InitialConditions`] set from the values at the position `index` of the arrays.
    fn initial_from_index(&self, index: usize) -> InitialConditions {
        InitialConditions {
            t0: 0.0,
            theta0: self.thetas[[index]],
            psip0: self.psips[[index]],
            rho0: self.rhos[[index]],
            zeta0: self.zetas[[index]],
            mu: self.mus[[index]],
        }
    }

    /// Creates a [`Particle`] set from the values at the position `index` of the arrays.
    fn particle_from_index(&self, index: usize) -> Particle {
        Particle::new(&self.initial_from_index(index))
    }

    /// Returns the length of the stored arrays.
    pub fn len(&self) -> usize {
        self.thetas.len()
    }

    /// Returns `true` if the arrays have a length of 0.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

array1D_getter_impl!(PoincareInit, thetas, thetas);
array1D_getter_impl!(PoincareInit, psips, psips);
array1D_getter_impl!(PoincareInit, rhos, rhos);
array1D_getter_impl!(PoincareInit, zetas, zetas);
array1D_getter_impl!(PoincareInit, mus, mus);

impl std::fmt::Debug for PoincareInit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PoincareInit")
            .field("Length", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_poincare_init_creation() {
        let p = PoincareInit::build(
            &[0.0, 1.0],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        assert_eq!(p.len(), 2);
        assert!(!p.is_empty());

        PoincareInit::build(
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap_err();
    }

    #[test]
    fn test_poincare_init_data_extraction() {
        let p = PoincareInit::build(
            &[0.0, 1.0],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        assert_eq!(p.thetas().len(), 2);
        assert_eq!(p.psips().len(), 2);
        assert_eq!(p.rhos().len(), 2);
        assert_eq!(p.zetas().len(), 2);
        assert_eq!(p.mus().len(), 2);
    }

    #[test]
    fn test_poincare_init_to_particles() {
        let p = PoincareInit::build(
            &[0.0, 1.0],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        let particles: Vec<Particle> = p.to_particles();
        assert_eq!(particles.len(), p.len());
    }
    #[test]
    fn test_poincare_init_dbg() {
        let p = PoincareInit::build(
            &[0.0, 1.0],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        let _ = format!("{:?}", p);
    }
}
