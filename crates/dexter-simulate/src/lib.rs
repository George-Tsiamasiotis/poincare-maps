//! Simulations of charged particles inside an [`Equilibrium`](dexter_equilibrium).
//!
//! ### [`Particle`]: A charged particle inside an [`Equilibrium`](dexter_equilibrium)
//!
//! The Particle corresponds to a proton with `m=q=1`. All its related quantities and calculated
//! time series are in *Normalized Units*.
//!
//! Routines:
//!
//! + [`Particle::integrate`]: Integrates the particle in a specific time interval.
//! + [`Particle::intersect`]: Integrates the particle, calculating its intersections with a
//!   constant `θ` or `ζ` surface.
//! + [`Particle::close`]: Integrates the particle, for a given number of `θ-ψ` periods, and calculates
//!   its `ωθ`, `ωζ` and `qkinetic`.
//! + [`Particle::classify`]: Classifies the particle's [`orbit`](crate::OrbitType) using its position
//!   on the `(E, Pζ, μ=const)` plane without integrating.
//!
//! ### [`Queue`]: A container for batch initializing and running Particle routines.
//!
//! A Queue is conveniently initialized from 1D [`arrays`](ndarray::ArrayBase) of initial conditions with
//! the help of the [`QueueInitialConditions`] helper struct.
//!
//! Routines:
//!
//! + [`Queue::integrate`]: Integrates the contained particles in a specific time interval.
//! + [`Queue::intersect`]: Integrates the contained particles, calculating their intersections with a
//!   constant `θ` or `ζ` surface. Otherwise known as a `Poincare map`.
//! + [`Queue::close`]: Integrates all the contained particles for a specific number of periods and
//!   calculates their frequencies.
//! + [`Queue::classify`]: Classifies all the contained particles' [`orbits`](crate::OrbitType), using
//!   their position on the (E, Pζ, μ=const) plane without integrating.
//! + [`Queue::classify_common_mu`]: An optimization to [`Queue::classify`] for classifying
//!   particles with common `μ`. Results to about 5-8 times better performance.
//!
//! ### COMs space
//!
//! The container type [`COMs`] provides methods and constructors for calculations on the `(E, Pζ, μ)`
//! space:
//!
//! + [`EnergyPzetaPlane`]: Representation of the COM space `(E, Pζ, μ=const)`.
//! + [`TrappedPassingBoundary`]: Representation of the Trapped-Passing boundary curves on the
//!   `(E, Pζ, μ=const)` space.
//! + [`energy_of_psi_grid`](COMs::energy_of_psi_grid): Calculation of the unperturbed Hamiltonian's
//!   value in a 2x2 `θ-ψ` grid
//! + [`energy_of_psip_grid`](COMs::energy_of_psip_grid): Calculation of the unperturbed Hamiltonian's
//!   value in a 2x2 `θ-ψp` grid
//!
//! ### Parallelism
//!
//! Using the [`rayon`] crate, particle routines can run in parallel. The number of threads to be
//! used can be adjusted with the [`set_num_threads`] function.
//!
//! ### Integration Configuration
//!
//! Each integration routine's set of parameters can be adjusted with the following structs:
//!
//! + [`SolverParams`]: Parameters passed to the solver.
//!
//! The integration is done with the RKF4(5) method. The step size can be adaptive by minimizing
//! the energy difference from step to step or minimizing the local truncation error (classic
//! RKF4(5)), or simple set to be constant.
//!
//! ### Constants of motion
//!
//! Calculations involving the Constants of Motion.
//!
//! + [`COMs`]: A container struct for the constants of motion in an unperturbed equilibrium.
//!   + [`COMs::energy_of_psi_grid`]: Calculation of the Energy in a 2D `θ-ψ` grid
//!   + [`COMs::energy_of_psip_grid`]: Calculation of the Energy in a 2D `θ-ψp` grid
//!

mod coms;
mod error;
mod particle;
mod queue;
mod solve;
mod state;

// ============== Public API

pub use dexter_common::{get_max_threads, set_num_threads};

pub use error::{COMError, SimulationError};

pub use solve::{FluxCoordinate, SolverParams, SteppingMethod};

pub use particle::{
    CoordinateSet, EnergyPzetaPosition, Frequencies, InitialConditions, InitialFlux,
    IntegrationStatus, IntersectParams, Intersection, OrbitType, Particle,
};

pub use queue::{Queue, QueueInitialConditions, Routine};
pub use queue::{poloidal_fluxes, toroidal_fluxes};

pub use coms::{COMs, EnergyPzetaPlane, TrappedPassingBoundary};

// ============== Configuration constants

/// Crate configuration constants.
pub mod constants {
    use super::SteppingMethod;

    // ================ Solver Parameters ================

    /// The default optimal step calculation method.
    pub const DEFAULT_STEPPING_METHOD: SteppingMethod = SteppingMethod::EnergyAdaptiveStep;

    /// The default maximum amount of steps a particle can make before terminating its integration.
    pub const DEFAULT_MAX_STEPS: usize = 1_000_000;

    /// The default initial time step for the RKF45 adaptive step method. The value is empirical.
    pub const DEFAULT_FIRST_STEP: f64 = 1e-1;

    /// The default safety factor of the solver. Should be less than 1.0.
    pub const DEFAULT_SAFETY_FACTOR: f64 = 0.9;

    /// The default relative tolerance of the energy difference in every step.
    pub const DEFAULT_ENERGY_REL_TOL: f64 = 1e-12;

    /// The default absolute tolerance of the energy difference in every step.
    pub const DEFAULT_ENERGY_ABS_TOL: f64 = 1e-14;

    /// The default relative tolerance of the local truncation error in every step.
    pub const DEFAULT_ERROR_REL_TOL: f64 = 1e-12;

    /// The default absolute tolerance of the local truncation error in every step.
    pub const DEFAULT_ERROR_ABS_TOL: f64 = 1e-14;

    // ============== COMs Space Parameters ==============

    /// The density of the trapped-passing boundary curves' points.
    ///
    /// A higher number is needed to better classify Potato and Stagnated orbits.
    pub const TRAPPED_PASSING_BOUNDARY_DENSITY: usize = 500;

    // ============== Poincare intersection ==============

    /// The maximum allowed relative difference between two consecutive poincare intersections.
    ///
    /// The difference must be smaller that the solver's truncation error.
    pub const ANGLE_INTERSECTION_THRESHOLD: f64 = 1e-9;

    // ================== Orbit closing ==================

    /// The default relative tolerance of the short-circuit `flux-flux0` check in the
    /// [`Particle::close`][crate::Particle::close] routine.
    ///
    /// This only affects performance, but should be smaller than [`FLUX_REL_TOL`].
    pub const SHORT_CIRCUIT_FLUX_REL_TOL: f64 = 1e-4;

    /// The default relative tolerance of the `flux-flux0` check in the
    /// [`Particle::close`][crate::Particle::close] routine.
    pub const FLUX_REL_TOL: f64 = 1e-6;
}
