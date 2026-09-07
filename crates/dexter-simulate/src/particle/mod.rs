//! Representation of a charged particle.

mod classify;
mod close;
mod evolution;
mod initial;
mod integrate;
mod intersect;

pub use classify::EnergyPzetaPosition;
pub use initial::{CoordinateSet, InitialConditions};
pub use intersect::{IntersectParams, Intersection};

// ===============================================================================================

use ndarray::Array1;
use rsl_interpolation::{Accelerator, Accelerator2d};
use std::time::Duration;

use dexter_common::export_array1D_getter_impl;
use dexter_machine::{DynModeCaches, Machine};

use crate::SolverParams;
use crate::coms::EnergyPzetaPlane;
use evolution::Evolution;

// ===============================================================================================

/// Container for the caching objects needed for the evaluations.
#[derive(Default, Debug)]
#[non_exhaustive]
pub(crate) struct IntegrationCaches {
    /// The 2D integration Accelerator.
    acc: Accelerator2d,
    /// The caches of the perturbation's modes.
    mode_caches: DynModeCaches,
}

impl IntegrationCaches {
    /// Returns a mutable reference to the 2D Accelerator.
    pub(crate) fn acc(&mut self) -> &mut Accelerator2d {
        &mut self.acc
    }

    /// Returns a mutable reference to the magnetic flux Accelerator.
    pub(crate) fn flux_acc(&mut self) -> &mut Accelerator {
        self.acc.xacc()
    }

    /// Returns a mutable reference to mode caches.
    pub(crate) fn mode_caches(&mut self) -> &mut DynModeCaches {
        &mut self.mode_caches
    }
}

// ===============================================================================================

/// A [`Particle`]'s integration status.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntegrationStatus {
    /// Initialized by [`InitialConditions`], not integrated.
    Initialized,
    /// [`InitialConditions`] have not been fully calculated yet.
    PartlyInitialized,
    /// Invalid [`InitialConditions`]. May occur when using Mixed variables with objects that
    /// cannot define them.
    InvalidInitialConditions,
    /// [`InitialConditions`] were out of bounds.
    OutOfBoundsInitialization,
    /// Reached the end of the integration successfully.
    Integrated,
    /// Intersections calculation successful.
    Intersected,
    /// Integrated for a certain amount of `θ-ψ` periods.
    ClosedPeriods(usize),
    /// Escaped the last closed flux surface (LCFS).
    Escaped,
    /// Escaped when performing a step on the modified system.
    ///
    /// This indicates that something is wrong in Hénon's trick implementation.
    ModStateEscaped,
    /// Calculated some intersections correctly but also timed out.
    IntersectedTimedOut,
    /// Calculated invalid intersections.
    InvalidIntersections,
    /// Timed out after a maximum number of steps.
    TimedOut(Duration),
    /// Simulation failed for unknown reasons.
    Failed(Box<str>),
}

// ===============================================================================================

/// A particle's calculated frequencies and `qkinetic`.
///
/// The frequencies are calculated through the [`Particle::close`] routine.
#[derive(Default, Clone)]
pub struct Frequencies {
    /// The particle's calculated `ωθ`.
    pub omega_theta: Option<f64>,
    /// The particle's calculated `ωζ`.
    pub omega_zeta: Option<f64>,
    /// The particle's calculated `qkinetic`.
    pub qkinetic: Option<f64>,
}

/// A particle's orbit type.
///
/// As described by [`R. B. White`], an orbit is classified depending on its location relative
/// to the well-defined `(E, Pζ, μ=const)`.
///
/// [`R. B. White`]: https://doi.org/10.1142/P440
#[derive(Default, Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum OrbitType {
    /// Particle has not been classified.
    #[default]
    Undefined,
    /// A Trapped-Lost particle.
    ///
    /// # Definition
    ///
    /// A particle is called *trapped* if there exists a mirror point where `ρ=0`.
    TrappedLost,
    /// A Trapped-Confined particle.
    ///
    /// # Definition
    ///
    /// A particle is called *trapped* if there exists a mirror point where `ρ=0`.
    TrappedConfined,
    /// A CoPassing-Lost particle.
    ///
    /// # Definition
    ///
    /// A particle is called *co-passing* if it is not trapped and it holds that `dot(θ)>0`.
    CoPassingLost,
    /// A CoPassing-Confined particle.
    ///
    /// # Definition
    ///
    /// A particle is called *co-passing* if it is not trapped and it holds that `dot(θ)>0`.
    CoPassingConfined,
    /// A CounterPassing-Lost particle.
    ///
    /// # Definition
    ///
    /// A particle is called *counter-passing* if it is not trapped and it holds that `dot(θ)<0`.
    CuPassingLost,
    /// A CounterPassing-Confined particle.
    ///
    /// # Definition
    ///
    /// A particle is called *counter-passing* if it is not trapped and it holds that `dot(θ)<0`.
    CuPassingConfined,
    /// A Potato particle.
    ///
    /// # Definition
    ///
    /// A particle's orbit is called a *potato* orbit if it is trapped but still circles the
    /// magnetic axis due to its drift. In the `(E, Pζ)` plane, those lie inside the intersection
    /// of the trapped-passing boundary and the magnetic axis parabola.
    Potato,
    /// A Potato particle.
    ///
    /// # Definition
    ///
    /// A particle is called *stagnated* if it always has positive parallel velocity but does
    /// not circle the magnetic axis. In the `(E, Pζ)` plane, those lie to the right of the
    /// trapped-passing boundary and above the magnetic axis parabola.
    Stagnated,
    /// Not falling under any of the other categories.
    Unclassified,
}

// ===============================================================================================
// ===============================================================================================

/// Representation of a charged particle.
pub struct Particle {
    /// The [`InitialConditions`] set of the particle.
    initial_conditions: InitialConditions,
    /// Status of the particle's integration.
    integration_status: IntegrationStatus,
    /// The time evolution of the particle.
    evolution: Evolution,
    /// The integration Accelerator and Mode caches.
    caches: IntegrationCaches,
    /// The particle's position on the `E-Pζ` plane, relative to the orbit classification curves.
    energy_pzeta_position: EnergyPzetaPosition,
    /// The particle's orbit type.
    orbit_type: OrbitType,
    /// The particle's calculated `ωθ`, `ωζ` and `qkinetic`.
    frequencies: Frequencies,
    /// The particle's initial Energy in Normalized Units. The Energy depends both on the initial
    /// conditions and the machine.
    initial_energy: Option<f64>,
    /// The particle's energy after the integration.
    final_energy: Option<f64>,
}

impl Particle {
    /// Creates a new [`Particle`] from a set of [`InitialConditions`].
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_simulate::*;
    /// let psip0 = InitialFlux::Poloidal(0.05);
    /// let initial = InitialConditions::boozer(0.0, psip0, 0.0, 3.14, 1e-5, 1e-6);
    /// let mut particle = Particle::new(&initial);
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    #[must_use]
    pub fn new(initial_conditions: &InitialConditions) -> Self {
        use CoordinateSet::*;
        let integration_status = match initial_conditions.coordinate_set() {
            BoozerToroidal | BoozerPoloidal => IntegrationStatus::Initialized, // Always succeeds
            MixedToroidal | MixedPoloidal => IntegrationStatus::PartlyInitialized, // Needs `finalize()`
        };
        Self {
            initial_conditions: initial_conditions.to_owned(),
            integration_status,
            evolution: Evolution::default(),
            caches: IntegrationCaches::default(),
            energy_pzeta_position: EnergyPzetaPosition::default(),
            orbit_type: OrbitType::default(),
            frequencies: Frequencies::default(),
            initial_energy: None,
            final_energy: None,
        }
    }
}

// Routines
impl Particle {
    /// Integrates the particle for a certain time interval.
    ///
    /// The time interval is in Normalized Units (inverse gyro-frequency).
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_simulate::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.05);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_perturbation(&perturbation)
    ///     .build();
    ///
    /// let psi0 = InitialFlux::Toroidal(0.015);
    /// let initial = InitialConditions::boozer(0.0, psi0, 0.0, 3.14, 1e-5, 1e-6);
    /// let mut particle = Particle::new(&initial);
    /// particle.integrate(machine, (0.0, 1e2), &SolverParams::default());
    ///
    /// assert_eq!(particle.integration_status(), IntegrationStatus::Integrated);
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    pub fn integrate(&mut self, machine: Machine, teval: (f64, f64), solver_params: &SolverParams) {
        integrate::integrate(self, machine, teval, solver_params);
    }

    /// Integrates the particle, calculating its intersections with a constant `θ` or `ζ` surface.
    ///
    /// The intersection surface, angle, and number of turns are configured with the helper struct
    /// [`IntersectParams`].
    ///
    /// Using the method described by [`Hénon`] we can force the solver to step exactly on the
    /// intersection surface.
    ///
    /// The differences between two consecutive values of the corresponding angle variable are
    /// guaranteed to be `2π +- ε`, where ε a number smaller than the solver's relative tolerance.
    ///
    /// [`Hénon`]: https://www.sciencedirect.com/science/article/abs/pii/0167278982900343
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_simulate::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.1);
    /// let qfactor = UnityQfactor::new(lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_perturbation(&perturbation)
    ///     .build();
    ///
    /// let psi0 = InitialFlux::Toroidal(0.02);
    /// let initial = InitialConditions::boozer(0.0, psi0, 3.14, 0.0, 1e-4, 1e-6);
    /// let mut particle = Particle::new(&initial);
    ///
    /// let intersect_params = IntersectParams::new(Intersection::ConstTheta, 3.14, 10);
    ///
    /// particle.intersect(machine, &intersect_params, &SolverParams::default());
    ///
    /// assert_eq!(particle.steps_stored(), 10);
    /// assert_eq!(particle.integration_status(), IntegrationStatus::Intersected);
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    pub fn intersect(
        &mut self,
        machine: Machine,
        intersect_params: &IntersectParams,
        solver_params: &SolverParams,
    ) {
        intersect::intersect(self, machine, intersect_params, solver_params);
    }

    /// Integrates the particle, for `periods` number of `θ-ψ` periods.
    ///
    /// When the particle approximately reaches its initial point after `periods` periods, it
    /// halts its integration and performs one more step using [`Hénon`]'s trick, to land
    /// itself on the initial point exactly, similar to how [`Particle::intersect`] lands exactly
    /// on the intersection surface.
    ///
    /// This routine also yields the particle's `ωθ`, `ωζ` and `qkinetic`.
    ///
    /// [`Hénon`]: https://www.sciencedirect.com/science/article/abs/pii/0167278982900343
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_simulate::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.05);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_perturbation(&perturbation)
    ///     .build();
    ///
    /// let psi0 = InitialFlux::Toroidal(0.015);
    /// let initial = InitialConditions::boozer(0.0, psi0, 0.0, 3.14, 1e-5, 1e-6);
    /// let mut particle = Particle::new(&initial);
    ///
    /// particle.close(machine, 1, &SolverParams::default());
    ///
    /// assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    pub fn close(&mut self, machine: Machine, periods: usize, solver_params: &SolverParams) {
        close::close(self, machine, periods, solver_params);
    }

    /// Classifies the particle's orbit using its position on the `(E, Pζ, μ=const)` plane without integrating.
    ///
    /// # Note
    ///
    /// This method is experimental. It is exact for LAR equilibria and approximately correct for
    /// tokamak equilibria, depending on how shaped they are.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_simulate::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.03);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
    /// ]);
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    ///
    /// let psi0 = InitialFlux::Toroidal(0.001);
    /// let pzeta0 = -0.8 * machine.qfactor().psip_last();
    /// let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);
    ///
    /// let mut particle = Particle::new(&initial);
    /// particle.classify(machine);
    ///
    /// assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    pub fn classify(&mut self, machine: Machine) {
        self._classify(machine, None);
    }

    /// Does the actual classification.
    ///
    /// OPTIM: Since the most common scenario is to classify particles with the same `μ`, we can
    /// generate the [`EnergyPzetaPlane`] only once and use it for all particles. This method should
    /// only be called by [`crate::Queue::classify_common_mu`].
    pub(crate) fn _classify(&mut self, objects: Machine, _plane: Option<&EnergyPzetaPlane>) {
        classify::classify(self, objects, _plane)
    }
}

// Getters
impl Particle {
    /// Returns the Particle's [`InitialConditions`].
    #[must_use]
    pub fn initial_conditions(&self) -> InitialConditions {
        self.initial_conditions.clone()
    }

    /// Returns the Particle's [`IntegrationStatus`].
    #[must_use]
    pub fn integration_status(&self) -> IntegrationStatus {
        self.integration_status.clone()
    }

    /// Returns the total number of steps taken during the integration.
    ///
    /// This number is not necessarily the same as the number of steps stored.
    #[must_use]
    pub fn steps_taken(&self) -> usize {
        self.evolution.steps_taken
    }

    /// Returns the number of steps stored in the time series arrays.
    #[must_use]
    pub fn steps_stored(&self) -> usize {
        self.evolution.steps_stored()
    }

    /// Returns the duration of the integration routine.
    #[must_use]
    pub fn duration(&self) -> Duration {
        self.evolution.duration
    }

    /// Returns the particle's initial energy.
    #[must_use]
    pub fn initial_energy(&self) -> Option<f64> {
        self.initial_energy
    }

    /// Returns the particle's final energy.
    #[must_use]
    pub fn final_energy(&self) -> Option<f64> {
        self.final_energy
    }

    /// Returns the variance of the energy array.
    #[must_use]
    pub fn energy_var(&self) -> Option<f64> {
        self.evolution.energy_var()
    }

    /// Returns the particle's [`EnergyPzetaPosition`].
    #[must_use]
    pub fn energy_pzeta_position(&self) -> EnergyPzetaPosition {
        self.energy_pzeta_position
    }

    /// Returns the particle's [`OrbitType`].
    #[must_use]
    pub fn orbit_type(&self) -> OrbitType {
        self.orbit_type.clone()
    }

    /// Returns the particle's [`Frequencies`].
    #[must_use]
    pub fn frequencies(&self) -> Frequencies {
        self.frequencies.clone()
    }

    /// Returns the particle's calculated `ωθ`.
    #[must_use]
    pub fn omega_theta(&self) -> Option<f64> {
        self.frequencies.omega_theta
    }

    /// Returns the particle's calculated `ωζ`.
    #[must_use]
    pub fn omega_zeta(&self) -> Option<f64> {
        self.frequencies.omega_zeta
    }

    /// Returns the particle's calculated `qkinetic`.
    #[must_use]
    pub fn qkinetic(&self) -> Option<f64> {
        self.frequencies.qkinetic
    }

    /// Prints the Accelerators' and mode caches' hits and misses.
    pub fn print_caches(&self) {
        println!("Particle caches {{");
        println!("\tFlux cache hits: {}", self.flux_cache_hits());
        println!("\tFlux cache misses: {}", self.flux_cache_misses());
        println!("\tTheta cache hits: {}", self.theta_cache_hits());
        println!("\tTheta cache misses: {}", self.theta_cache_misses());
        println!("\tMode cache hits: {}", self.mode_cache_hits());
        println!("\tMode cache misses: {}", self.mode_cache_misses());
        println!("}}");
    }

    /// Discares the time evolution arrays, keeping metadata such as duration, step count, etc.
    pub fn discard_vecs(&mut self) {
        self.evolution.discard_vecs();
    }

    /// Stores the integration caches in the particle. Should be called after every integration routine.
    pub(crate) fn store_caches(&mut self, caches: IntegrationCaches) {
        self.caches = caches
    }

    /// Returns the Accelerator's magnetic flux cache hits.
    #[must_use]
    pub fn flux_cache_hits(&self) -> usize {
        self.caches.acc.clone().xacc().hits()
    }

    /// Returns the Accelerator's magnetic flux cache misses.
    #[must_use]
    pub fn flux_cache_misses(&self) -> usize {
        self.caches.acc.clone().xacc().misses()
    }

    /// Returns the Accelerator's theta angle cache hits.
    #[must_use]
    pub fn theta_cache_hits(&self) -> usize {
        self.caches.acc.clone().yacc().hits()
    }

    /// Returns the Accelerator's theta angle cache misses.
    #[must_use]
    pub fn theta_cache_misses(&self) -> usize {
        self.caches.acc.clone().yacc().misses()
    }

    /// Returns the mode cache's hits.
    #[must_use]
    pub fn mode_cache_hits(&self) -> usize {
        self.caches
            .mode_caches
            .iter()
            .fold(0, |acc, cache| acc + cache.hits())
    }

    /// Returns the mode cache's misses.
    #[must_use]
    pub fn mode_cache_misses(&self) -> usize {
        self.caches
            .mode_caches
            .iter()
            .fold(0, |acc, cache| acc + cache.misses())
    }

    export_array1D_getter_impl!(t_array, evolution, t);
    export_array1D_getter_impl!(psi_array, evolution, psi_array);
    export_array1D_getter_impl!(psip_array, evolution, psip_array);
    export_array1D_getter_impl!(theta_array, evolution, theta_array);
    export_array1D_getter_impl!(zeta_array, evolution, zeta_array);
    export_array1D_getter_impl!(rho_array, evolution, rho_array);
    export_array1D_getter_impl!(mu_array, evolution, mu_array);
    export_array1D_getter_impl!(ptheta_array, evolution, ptheta_array);
    export_array1D_getter_impl!(pzeta_array, evolution, pzeta_array);
    export_array1D_getter_impl!(energy_array, evolution, energy_array);
}

impl Clone for Particle {
    fn clone(&self) -> Self {
        Self {
            initial_conditions: self.initial_conditions.clone(),
            integration_status: self.integration_status.clone(),
            evolution: self.evolution.clone(),
            caches: IntegrationCaches::default(),
            energy_pzeta_position: self.energy_pzeta_position,
            orbit_type: self.orbit_type.clone(),
            frequencies: self.frequencies.clone(),
            initial_energy: self.initial_energy,
            final_energy: self.final_energy,
        }
    }
}

impl std::fmt::Debug for Frequencies {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fn stringify(field: Option<f64>) -> String {
            if let Some(value) = field {
                format!("{value:.7}")
            } else {
                String::from("Not calculated")
            }
        }

        f.debug_struct("Frequencies")
            .field("omega_theta", &stringify(self.omega_theta))
            .field("omega_zeta", &stringify(self.omega_zeta))
            .field("qkinetic", &stringify(self.qkinetic))
            .finish()
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("initial conditions", &self.initial_conditions)
            .field("integration status", &self.integration_status)
            .field("evolution", &self.evolution)
            .field("orbit_type", &self.orbit_type)
            .field("frequencies", &self.frequencies)
            .field("initial energy", &self.initial_energy.unwrap_or(f64::NAN))
            .field("final energy  ", &self.final_energy.unwrap_or(f64::NAN))
            .field("energy variance", &self.energy_var().unwrap_or(f64::NAN))
            .finish()
    }
}
