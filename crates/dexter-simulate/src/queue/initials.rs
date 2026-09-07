//! Queue Initial Conditions for batch Particle initialization.

use dexter_machine::MagneticFlux;
use ndarray::{Array1, ArrayView1};

use crate::{CoordinateSet, InitialConditions, Particle, SimulationError};

/// Sets of initial conditions for initializing a [`Queue`](crate::Queue`).
///
/// The [`toroidal_fluxes`] and [`poloidal_fluxes`] helper functions can be used to create the
/// initial fluxes array.
#[non_exhaustive]
#[derive(Clone)]
pub struct QueueInitialConditions {
    /// The initial times `t0`.
    t0: Array1<f64>,
    /// The initial `fluxes`.
    flux0: Array1<MagneticFlux>,
    /// The initial `thetas`.
    theta0: Array1<f64>,
    /// The initial `zetas`.
    zeta0: Array1<f64>,
    /// The initial `mus`.
    mu0: Array1<f64>,
    /// The initial `rhos`.
    rho0: Option<Array1<f64>>,
    /// The initial `Pzetas`.
    pzeta0: Option<Array1<f64>>,
    /// Which coordinate set were the initial conditions defined.
    coordinate_set: CoordinateSet,
}

impl QueueInitialConditions {
    /// Creates a new [`QueueInitialConditions`] in Boozer coordinates.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// use MagneticFlux::*;
    /// let initial_conditions = QueueInitialConditions::boozer(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[1e-3, 2e-3],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// # Ok::<_, SimulationError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the inputs do not have the same length or are empty.
    pub fn boozer(
        t0: &[f64],
        flux0: &[MagneticFlux],
        theta0: &[f64],
        zeta0: &[f64],
        rho0: &[f64],
        mu0: &[f64],
    ) -> Result<Self, SimulationError> {
        let len = t0.len();
        if !(flux0.len() == len
            && theta0.len() == len
            && zeta0.len() == len
            && rho0.len() == len
            && mu0.len() == len)
        {
            return Err(SimulationError::QueueInitialConditionsMismatch);
        }
        let coordinate_set = match flux0
            .first()
            .copied()
            .ok_or(SimulationError::QueueInitialConditionsEmptyInput)?
        {
            MagneticFlux::Toroidal(_) => CoordinateSet::BoozerToroidal,
            MagneticFlux::Poloidal(_) => CoordinateSet::BoozerPoloidal,
        };
        Ok(Self {
            t0: Array1::from(t0.to_vec()),
            flux0: Array1::from(flux0.to_vec()),
            theta0: Array1::from(theta0.to_vec()),
            zeta0: Array1::from(zeta0.to_vec()),
            mu0: Array1::from(mu0.to_vec()),
            rho0: Some(Array1::from(rho0.to_vec())),
            pzeta0: None,
            coordinate_set,
        })
    }

    /// Creates a new [`QueueInitialConditions`] in Mixed coordinates.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// use MagneticFlux::*;
    /// let initial_conditions = QueueInitialConditions::mixed(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[-0.025, 0.027],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// # Ok::<_, SimulationError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the inputs do not have the same length or are empty.
    pub fn mixed(
        t0: &[f64],
        flux0: &[MagneticFlux],
        theta0: &[f64],
        zeta0: &[f64],
        pzeta0: &[f64],
        mu0: &[f64],
    ) -> Result<Self, SimulationError> {
        let len = t0.len();
        if !(flux0.len() == len
            && theta0.len() == len
            && zeta0.len() == len
            && pzeta0.len() == len
            && mu0.len() == len)
        {
            return Err(SimulationError::QueueInitialConditionsMismatch);
        }
        let coordinate_set = match flux0
            .first()
            .copied()
            .ok_or(SimulationError::QueueInitialConditionsEmptyInput)?
        {
            MagneticFlux::Toroidal(_) => CoordinateSet::MixedToroidal,
            MagneticFlux::Poloidal(_) => CoordinateSet::MixedPoloidal,
        };
        Ok(Self {
            t0: Array1::from(t0.to_vec()),
            flux0: Array1::from(flux0.to_vec()),
            theta0: Array1::from(theta0.to_vec()),
            zeta0: Array1::from(zeta0.to_vec()),
            mu0: Array1::from(mu0.to_vec()),
            rho0: None,
            pzeta0: Some(Array1::from(pzeta0.to_vec())),
            coordinate_set,
        })
    }

    /// Creates a new `QueueInitialConditions` from a slice of particles by extracting their
    /// `InitialConditions`.
    ///
    /// The particles do not have to be fully initialized.
    ///
    /// # Panics
    ///
    /// Panics if *not* all the particles are defined in the same coordinate set (boozer/mixed).
    pub(crate) fn from_particles(particles: &[Particle]) -> Self {
        let mut t0 = Vec::with_capacity(particles.len());
        let mut flux0 = Vec::with_capacity(particles.len());
        let mut theta0 = Vec::with_capacity(particles.len());
        let mut zeta0 = Vec::with_capacity(particles.len());
        let mut mu0 = Vec::with_capacity(particles.len());

        // If any values are found, store them here, and the create the array
        let mut rho0 = Vec::new();
        let mut pzeta0 = Vec::new();

        for particle in particles {
            let initial = particle.initial_conditions();
            t0.push(initial.t0);
            flux0.push(initial.flux0);
            theta0.push(initial.theta0);
            zeta0.push(initial.zeta0);
            mu0.push(initial.mu0);
            match initial.coordinate_set() {
                CoordinateSet::BoozerToroidal | CoordinateSet::BoozerPoloidal => {
                    rho0.push(initial.rho0.expect("is defined"))
                }
                CoordinateSet::MixedToroidal | CoordinateSet::MixedPoloidal => {
                    pzeta0.push(initial.pzeta0.expect("is defined"))
                }
            }
        }

        // Sanity check
        let len = t0.len();
        let other_len = rho0.len().max(pzeta0.len());
        if !(flux0.len() == len
            && theta0.len() == len
            && zeta0.len() == len
            && mu0.len() == len
            && other_len == len)
        {
            unreachable!()
        }

        Self {
            t0: Array1::from(t0),
            flux0: Array1::from(flux0),
            theta0: Array1::from(theta0),
            zeta0: Array1::from(zeta0),
            mu0: Array1::from(mu0),
            rho0: (!rho0.is_empty()).then(|| Array1::from(rho0)),
            pzeta0: (!pzeta0.is_empty()).then(|| Array1::from(pzeta0)),
            coordinate_set: particles[0].initial_conditions().coordinate_set(),
        }
    }

    /// Creates a [`Vec<Particle>`] initialized from each set of the contained
    /// [`InitialConditions`].
    #[must_use]
    pub fn to_particles(&self) -> Vec<Particle> {
        let mut particles = Vec::with_capacity(self.len());
        for index in 0..self.len() {
            particles.push(self.particle_from_index(index));
        }
        particles
    }

    /// Creates a [`Particle`] by picking the [`InitialConditions`] at the `index` position of the arrays.
    #[must_use]
    pub fn particle_from_index(&self, index: usize) -> Particle {
        Particle::new(&self.initial_from_index(index))
    }

    /// Creates an [`InitialConditions`] from the `index` position of the arrays.
    #[expect(clippy::missing_panics_doc, reason = "safe")]
    #[must_use]
    pub fn initial_from_index(&self, index: usize) -> InitialConditions {
        use CoordinateSet::*;
        match self.coordinate_set {
            BoozerToroidal | BoozerPoloidal => InitialConditions::boozer(
                self.t0[[index]],
                self.flux0[[index]],
                self.theta0[[index]],
                self.zeta0[[index]],
                self.rho0.as_ref().expect("always `Some` in `boozer`")[[index]],
                self.mu0[[index]],
            ),
            MixedToroidal | MixedPoloidal => InitialConditions::mixed(
                self.t0[[index]],
                self.flux0[[index]],
                self.theta0[[index]],
                self.zeta0[[index]],
                self.pzeta0.as_ref().expect("always `Some` in `mixed`")[[index]],
                self.mu0[[index]],
            ),
        }
    }

    /// Returns the number of [`InitialConditions`] sets.
    #[must_use]
    pub fn len(&self) -> usize {
        self.t0.len()
    }

    /// Return `true` if there are no [`InitialConditions`] sets.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Helper function to create an array of `ψ0`s, to be passed to
/// [`QueueInitialConditions`]' constructors.
///
/// # Example
/// ```
/// # use dexter_simulate::*;
/// # use ndarray::Array1;
/// let psi0s = toroidal_fluxes(&Array1::linspace(0.0, 0.45, 100).to_vec());
/// ```
#[must_use]
pub fn toroidal_fluxes(values: &[f64]) -> Array1<MagneticFlux> {
    Array1::from_iter(values.iter().map(|value| MagneticFlux::Toroidal(*value)))
}

/// Helper function to create an array of `ψp0`s, to be passed to
/// [`QueueInitialConditions`]' constructors.
///
/// # Example
/// ```
/// # use dexter_simulate::*;
/// # use ndarray::Array1;
/// let psip0s = poloidal_fluxes(&Array1::linspace(0.0, 0.45, 100).to_vec());
/// ```
#[must_use]
pub fn poloidal_fluxes(values: &[f64]) -> Array1<MagneticFlux> {
    Array1::from_iter(values.iter().map(|value| MagneticFlux::Poloidal(*value)))
}

/// Getters.
impl QueueInitialConditions {
    /// Returns the initial times array.
    #[must_use]
    pub fn t_array(&self) -> Array1<f64> {
        self.t0.clone()
    }

    /// Returns the initial [`MagneticFlux`] array.
    #[must_use]
    pub fn flux_array(&self) -> Array1<MagneticFlux> {
        self.flux0.clone()
    }

    /// Returns the initial `θ` array.
    #[must_use]
    pub fn theta_array(&self) -> Array1<f64> {
        self.theta0.clone()
    }

    /// Returns the initial `ζ` array.
    #[must_use]
    pub fn zeta_array(&self) -> Array1<f64> {
        self.zeta0.clone()
    }

    /// Returns the initial `μ` array.
    #[must_use]
    pub fn mu_array(&self) -> Array1<f64> {
        self.mu0.clone()
    }

    /// Returns an [`ArrayView1`] of the `mu` array.
    #[must_use]
    pub(crate) fn mu_array_view(&self) -> ArrayView1<'_, f64> {
        self.mu0.view()
    }

    /// Returns the initial `ρ` array.
    #[must_use]
    pub fn rho_array(&self) -> Option<Array1<f64>> {
        self.rho0.clone()
    }

    /// Returns the initial `Pζ` array.
    #[must_use]
    pub fn pzeta_array(&self) -> Option<Array1<f64>> {
        self.pzeta0.clone()
    }
}

impl std::fmt::Debug for QueueInitialConditions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("QueueInitialConditions")
            .field("Length", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use crate::IntegrationStatus;

    use super::*;
    use MagneticFlux::*;

    #[test]
    fn queue_boozer_initial_conditions_creation() {
        let initial = QueueInitialConditions::boozer(
            &[0.0, 0.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[0.0, 1.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        assert_eq!(initial.len(), 2);
        assert!(!initial.is_empty());
        assert_eq!(initial.coordinate_set, CoordinateSet::BoozerToroidal);
        let _ = format!("{:?}", initial);

        let _ = QueueInitialConditions::boozer(
            &[0.0, 1.0, 2.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap_err();
    }

    #[test]
    fn queue_mixed_initial_conditions_creation() {
        let initial = QueueInitialConditions::mixed(
            &[0.0, 0.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[0.0, 1.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        assert_eq!(initial.len(), 2);
        assert!(!initial.is_empty());
        assert_eq!(initial.coordinate_set, CoordinateSet::MixedToroidal);
        let _ = format!("{:?}", initial);

        let _ = QueueInitialConditions::mixed(
            &[0.0, 1.0, 2.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap_err();
    }

    #[test]
    fn queue_boozer_initial_conditions_data_extraction() {
        let initial = QueueInitialConditions::boozer(
            &[0.0, 1.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        assert_eq!(initial.t_array().len(), 2);
        assert_eq!(initial.flux_array().len(), 2);
        assert_eq!(initial.theta_array().len(), 2);
        assert_eq!(initial.zeta_array().len(), 2);
        assert_eq!(initial.mu_array().len(), 2);
        assert_eq!(initial.rho_array().unwrap().len(), 2);
        assert!(initial.pzeta_array().is_none());
    }

    #[test]
    fn queue_mixed_initial_conditions_data_extraction() {
        let initial = QueueInitialConditions::mixed(
            &[0.0, 1.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        assert_eq!(initial.t_array().len(), 2);
        assert_eq!(initial.flux_array().len(), 2);
        assert_eq!(initial.theta_array().len(), 2);
        assert_eq!(initial.zeta_array().len(), 2);
        assert_eq!(initial.mu_array().len(), 2);
        assert_eq!(initial.pzeta_array().unwrap().len(), 2);
        assert!(initial.rho_array().is_none());
    }

    #[test]
    fn queue_boozer_initial_conditions_to_particles() {
        let initial = QueueInitialConditions::boozer(
            &[0.0, 1.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        let particles: Vec<Particle> = initial.to_particles();
        particles.iter().for_each(|particle| {
            assert_eq!(
                particle.integration_status(),
                IntegrationStatus::Initialized
            )
        });
        assert_eq!(particles.len(), initial.len());
    }

    #[test]
    fn queue_mixed_initial_conditions_to_particles() {
        let initial = QueueInitialConditions::mixed(
            &[0.0, 1.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        let particles: Vec<Particle> = initial.to_particles();
        particles.iter().for_each(|particle| {
            assert_eq!(
                particle.integration_status(),
                IntegrationStatus::PartlyInitialized
            )
        });
        assert_eq!(particles.len(), initial.len());
    }

    #[test]
    fn fluxes_queue_init_helpers() {
        let _: Array1<MagneticFlux> = toroidal_fluxes(&[0.1, 0.2]);
        let _: Array1<MagneticFlux> = poloidal_fluxes(&[0.1, 0.2]);
    }

    #[test]
    fn from_boozer_particles() {
        let psi0 = Toroidal(0.001);
        let initial1 = InitialConditions::boozer(0.0, psi0, 1.0, 0.0, 0.1, 0.0);
        let initial2 = InitialConditions::boozer(0.0, psi0, 2.0, 0.0, 0.1, 0.0);
        let particle1 = Particle::new(&initial1);
        let particle2 = Particle::new(&initial2);

        let initials = QueueInitialConditions::from_particles(&[particle1, particle2]);
        let queue = crate::Queue::new(&initials);

        assert_eq!(initial1.theta0, queue[0].initial_conditions().theta0);
        assert_eq!(initial2.theta0, queue[1].initial_conditions().theta0);
    }

    #[test]
    fn from_mixed_particles() {
        let psi0 = Toroidal(0.001);
        let initial1 = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, 0.0, 0.0);
        let initial2 = InitialConditions::mixed(0.0, psi0, 2.0, 0.0, 0.0, 0.0);
        let particle1 = Particle::new(&initial1);
        let particle2 = Particle::new(&initial2);

        let initials = QueueInitialConditions::from_particles(&[particle1, particle2]);
        let queue = crate::Queue::new(&initials);

        assert_eq!(initial1.theta0, queue[0].initial_conditions().theta0);
        assert_eq!(initial2.theta0, queue[1].initial_conditions().theta0);
    }

    #[test]
    #[should_panic]
    fn from_boozer_mixed_particles() {
        let psi0 = Toroidal(0.001);
        let initial1 = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, 0.0, 0.0);
        let initial2 = InitialConditions::boozer(0.0, psi0, 2.0, 0.0, 0.1, 0.0);
        let particle1 = Particle::new(&initial1);
        let particle2 = Particle::new(&initial2);

        let initials = QueueInitialConditions::from_particles(&[particle1, particle2]);
        let queue = crate::Queue::new(&initials);

        assert_eq!(initial1.theta0, queue[0].initial_conditions().theta0);
        assert_eq!(initial2.theta0, queue[1].initial_conditions().theta0);
    }
}
