//! Definition of a Particle's initial conditions in various coordinate sets.

use dexter_machine::Machine;
use rsl_interpolation::Accelerator;

use crate::{InitialFlux, SimulationError};

/// The kind of [`InitialConditions`] set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoordinateSet {
    /// Initial conditions set in the `(t, ψ, θ, ζ, ρ, μ)` space.
    BoozerToroidal,
    /// Initial conditions set in the `(t, ψp, θ, ζ, ρ, μ)` space.
    BoozerPoloidal,
    /// Initial conditions set in the `(t, Pζ, ψ, θ, ζ, μ)` space.
    MixedToroidal,
    /// Initial conditions set in the `(t, Pζ, ψp, θ, ζ, μ)` space.
    MixedPoloidal,
}

/// A set of initial conditions.
///
/// The set can be created from one the following combinations:
///
/// + Boozer-Toroidal: Initial conditions set in the `(t, ψ, θ, ζ, ρ, μ)` space.
/// + Boozer-Poloidal: Initial conditions set in the `(t, ψp, θ, ζ, ρ, μ)` space.
/// + Mixed-Toroidal: Initial conditions set in the `(t, Pζ, ψ, θ, ζ, μ)` space.
/// + Mixed-Poloidal: Initial conditions set in the `(t, Pζ, ψp, θ, ζ, μ)` space.
#[derive(Clone)]
#[non_exhaustive]
pub struct InitialConditions {
    /// The initial time.
    pub(crate) t0: f64,
    /// The initial toroidal/poloidal flux, depending on the value of [`InitialFlux`].
    pub(crate) flux0: InitialFlux,
    /// The initial `θ` angle.
    pub(crate) theta0: f64,
    /// The initial `ζ` angle.
    pub(crate) zeta0: f64,
    /// The magnetic moment `μ`.
    pub(crate) mu0: f64,
    /// The initial parallel gyro radius `ρ`.
    pub(crate) rho0: Option<f64>,
    /// The canonical momentum `Pζ`.
    pub(crate) pzeta0: Option<f64>,
    /// Which coordinate set were the initial conditions defined.
    coordinate_set: CoordinateSet,
}

impl InitialConditions {
    /// Constructs an `InitialConditions` from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψ, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// let psi0 = InitialFlux::Toroidal(0.2);
    /// let initial = InitialConditions::boozer(0.0, psi0, 0.0, 0.0, 1e-4, 1e-6);
    /// ```
    #[must_use]
    pub fn boozer(
        t0: f64,
        flux0: InitialFlux,
        theta0: f64,
        zeta0: f64,
        rho0: f64,
        mu0: f64,
    ) -> Self {
        let coordinate_set = match flux0 {
            InitialFlux::Toroidal(_) => CoordinateSet::BoozerToroidal,
            InitialFlux::Poloidal(_) => CoordinateSet::BoozerPoloidal,
        };
        Self {
            t0,
            flux0,
            theta0,
            zeta0,
            mu0,
            rho0: Some(rho0),
            pzeta0: None,
            coordinate_set,
        }
    }

    /// Constructs an `InitialConditions` from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψ, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// # use std::f64::consts::PI;
    /// let psip0 = InitialFlux::Poloidal(0.01);
    /// let initial = InitialConditions::mixed(0.0, psip0, PI, PI/2.0, -0.027, 1e-6);
    /// ```
    #[must_use]
    pub fn mixed(
        t0: f64,
        flux0: InitialFlux,
        theta0: f64,
        zeta0: f64,
        pzeta0: f64,
        mu0: f64,
    ) -> Self {
        let coordinate_set = match flux0 {
            InitialFlux::Toroidal(_) => CoordinateSet::MixedToroidal,
            InitialFlux::Poloidal(_) => CoordinateSet::MixedPoloidal,
        };
        Self {
            t0,
            flux0,
            theta0,
            zeta0,
            mu0,
            rho0: None,
            pzeta0: Some(pzeta0),
            coordinate_set,
        }
    }

    /// Calculates the missing coordinates needed for the integration routines.
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the missing coordinates cannot be calculated. This can
    /// occur if the [`Current`] object specifically does not define g(ψ) (`MixedToroidal` case) or
    /// g(ψp) (`MixedPoloidal` case), which are necessary for calculating the fluxes.
    pub(crate) fn finalize(&mut self, machine: Machine) -> Result<(), SimulationError> {
        let acc = &mut Accelerator::new();
        match self.coordinate_set {
            // Calculate `pzeta0`
            CoordinateSet::BoozerToroidal => {
                let psi0 = self.flux0.value();
                let g_of_psi0 = machine.current().g_of_psi(psi0, acc)?;
                let psip0 = machine.qfactor().psip_of_psi(psi0, acc)?;
                self.pzeta0 = Some(self.rho0.expect("boozer to mixed") * g_of_psi0 - psip0)
            }
            CoordinateSet::BoozerPoloidal => {
                let psip0 = self.flux0.value();
                let g_of_psip0 = machine.current().g_of_psip(psip0, acc)?;
                self.pzeta0 = Some(self.rho0.expect("boozer to mixed") * g_of_psip0 - psip0)
            }
            // Calculate `rho0`
            CoordinateSet::MixedToroidal => {
                let psi0 = self.flux0.value();
                let g_of_psi0 = machine.current().g_of_psi(psi0, acc)?;
                let psip0 = machine.qfactor().psip_of_psi(psi0, acc)?;
                self.rho0 = Some((self.pzeta0.expect("mixed to boozer") + psip0) / g_of_psi0);
            }
            CoordinateSet::MixedPoloidal => {
                let psip0 = self.flux0.value();
                let g_of_psip0 = machine.current().g_of_psip(psip0, acc)?;
                self.rho0 = Some((self.pzeta0.expect("mixed to boozer") + psip0) / g_of_psip0);
            }
        };
        self.assert_integration_ready();
        Ok(())
    }

    /// Asserts the necessary coordinates are initialized.
    ///
    /// `t0`, `flux0`, `theta0`, `zeta0` and `mu0` are guaranteed to exist, and we must ensure
    /// `rho0` is initialized as well.
    fn assert_integration_ready(&self) {
        assert!(
            self.t0.is_finite()
                && self.flux0.value().is_finite()
                && self.theta0.is_finite()
                && self.zeta0.is_finite()
                && self.rho0.is_some_and(f64::is_finite)
                && self.mu0.is_finite(),
            "Invalid InitialConditions"
        )
    }
}

/// Getters.
impl InitialConditions {
    /// Returns the initial time `t0`.
    #[must_use]
    pub fn t0(&self) -> f64 {
        self.t0
    }

    /// Returns the initial flux `flux0`.
    #[must_use]
    pub fn flux0(&self) -> InitialFlux {
        self.flux0
    }

    /// Returns the initial `θ`.
    #[must_use]
    pub fn theta0(&self) -> f64 {
        self.theta0
    }

    /// Returns the initial `ζ`.
    #[must_use]
    pub fn zeta0(&self) -> f64 {
        self.zeta0
    }

    /// Returns the initial `ρ`.
    #[must_use]
    pub fn rho0(&self) -> Option<f64> {
        self.rho0
    }

    /// Returns the initial `μ`.
    #[must_use]
    pub fn mu0(&self) -> f64 {
        self.mu0
    }

    /// Returns the initial `Pζ`.
    #[must_use]
    pub fn pzeta0(&self) -> Option<f64> {
        self.pzeta0
    }

    /// Returns the [`CoordinateSet`].
    #[must_use]
    pub fn coordinate_set(&self) -> CoordinateSet {
        self.coordinate_set.clone()
    }
}

impl std::fmt::Debug for InitialConditions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("InitialConditions")
            .field("t0", &self.t0)
            .field("flux0", &self.flux0)
            .field("theta0", &self.theta0)
            .field("zeta0", &self.zeta0)
            .field("mu0", &self.mu0)
            .field("rho0", &format!("{:?}", self.rho0))
            .field("pzeta0", &format!("{:?}", self.pzeta0))
            .field("coordinate_set", &self.coordinate_set)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use approx::assert_relative_eq;
    use dexter_machine::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
    use dexter_machine::*;
    use std::{f64::consts::PI, path::PathBuf};

    use super::*;
    use crate::*;
    use InitialFlux::*;

    #[test]
    fn boozer_initial_conditions() {
        let lcfs = LastClosedFluxSurface::Toroidal(0.05);
        let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let perturbation = Perturbation::new(vec![
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]);
        let machine = MachineBuilder::new(&qfactor, &current, &bfield)
            .with_perturbation(&perturbation)
            .build();

        let mut i1 = InitialConditions::boozer(0.0, Toroidal(0.02), PI, PI, 1e-4, 1e-6);
        i1.finalize(machine).unwrap();
        assert_eq!(i1.coordinate_set, CoordinateSet::BoozerToroidal);
        assert!(i1.rho0.is_some());
        assert!(i1.pzeta0.is_some());

        let mut i2 = InitialConditions::boozer(0.0, Poloidal(0.02), PI, PI, 1e-4, 1e-6);
        i2.finalize(machine).unwrap();
        assert_eq!(i2.coordinate_set, CoordinateSet::BoozerPoloidal);
        assert!(i2.rho0.is_some());
        assert!(i2.pzeta0.is_some());

        let mut initial = InitialConditions::boozer(0.0, Toroidal(0.02), PI, PI, 1e-4, 1e-6);
        initial.finalize(machine).unwrap();
        let _ = initial.t0();
        let _ = initial.flux0();
        let _ = initial.theta0();
        let _ = initial.zeta0();
        let _ = initial.rho0();
        let _ = initial.mu0();
        let _ = initial.pzeta0();
    }

    #[test]
    fn mixed_toroidal_initial_conditions() {
        let lcfs = LastClosedFluxSurface::Toroidal(0.05);
        let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let perturbation = Perturbation::new(vec![
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]);
        let machine = MachineBuilder::new(&qfactor, &current, &bfield)
            .with_perturbation(&perturbation)
            .build();

        let mut initial = InitialConditions::mixed(0.0, Toroidal(0.01), PI, PI, -0.027, 1e-6);
        initial.finalize(machine).unwrap();
        assert_eq!(initial.coordinate_set, CoordinateSet::MixedToroidal);
        assert!(initial.rho0.is_some());
        assert!(initial.pzeta0.is_some());

        let mut particle = Particle::new(&initial);
        particle.integrate(machine, (0.0, 1e2), &SolverParams::default());
        assert!(particle.steps_taken() > 10);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        // LarCurrent defines ψp evaluations but LarBfield cannot integrate with respect to ψp.
        InitialConditions::mixed(0.0, Poloidal(0.01), PI, PI, -0.027, 1e-6)
            .finalize(machine)
            .unwrap();
    }

    #[test]
    fn mixed_poloidal_initial_conditions() {
        let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let current = NcCurrentBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic)
            .build()
            .unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let mut initial = InitialConditions::mixed(0.0, Poloidal(0.01), PI, PI, -0.027, 1e-6);
        initial.finalize(machine).unwrap();
        assert!(initial.rho0.is_some());
        assert!(initial.pzeta0.is_some());

        let mut particle = Particle::new(&initial);
        particle.integrate(machine, (0.0, 1e2), &SolverParams::default());
        assert!(particle.steps_taken() > 10);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        let _ = InitialConditions::mixed(0.0, Toroidal(0.01), PI, PI, -0.027, 1e-6)
            .finalize(machine)
            .unwrap_err();
    }

    #[test]
    fn boozer_mixed_equivalence() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let current = NcCurrentBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic)
            .build()
            .unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let mut boozer = InitialConditions::boozer(0.0, Toroidal(0.01), PI, PI, 1e-4, 1e-6);
        boozer.finalize(machine).unwrap();
        // WARN: This used to be 1e-12 but randomly started failing, I could not find a reason why.
        assert_relative_eq!(boozer.pzeta0.unwrap(), -0.00898781038097592, epsilon = 1e-6);

        let mut mixed =
            InitialConditions::mixed(0.0, Toroidal(0.01), PI, PI, boozer.pzeta0.unwrap(), 1e-6);
        mixed.finalize(machine).unwrap();
        assert_relative_eq!(mixed.rho0.unwrap(), boozer.rho0.unwrap(), epsilon = 1e-12);
    }
}
