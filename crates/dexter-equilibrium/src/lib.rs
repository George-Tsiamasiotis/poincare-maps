#![doc = include_str!("../README.md")]
//!
//! # Equilibrium objects
//!
//! + Representations of an equilibrium's general geometry. Provides interpolation methods between `ψ`, `ψp`, `r`, `R`, `Z`, `J`.
//!     - [`LarGeometry`]: Analytical Large Aspect Ratio Geometry of a circular device.
//!     - [`NcGeometry`]: Geometry of the netCDF equilibrium
//!
//! + Representations of the q-factor profile:
//!     - [`UnityQfactor`]: q-factor profile of q = 1 and ψ=ψp.
//!     - [`ParabolicQfactor`]: q-factor of parabolic q(ψ) profile.
//!     - [`NcQfactor`]: q-factor reconstructed from a netCDF file.
//!
//! + Representations of the plasma currents:
//!     - [`LarCurrent`]: Large Aspect Ration plasma current with g=1 and I=0.
//!     - [`NcCurrent`]: Plasma current reconstructed from a netCDF file.
//!
//! + Representations of the magnetic field:
//!     - [`LarBfield`]: Large Aspect Ratio magnetic field with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).
//!     - [`NcBfield`]: Magnetic reconstructed from a netCDF file.
//!
//! + Representations of single perturbation modes:
//!     - [`FluteMode`]: Single analytical flute mode of the form `α*cos(mθ-nζ+φ)`.
//!     - [`NcFluteMode`]: Single numerical flute mode from a netCDF file of the form
//!     `α(ψ/ψp) * cos(mθ-nζ+φ(ψ/ψp))`
//!
//! + Representations of Perturbations.
//!     - [`Perturbation`]: A sum of an arbitrary number of [`Modes`](Mode).
//!
//! ## Evaluations:
//!
//! + [`Geometry`]: Conversions to laboratory quantities.
//! + [`FluxCommute`]: Conversion between the two flux coordinates `ψ` and `ψp`
//! + [`Qfactor`]: Evaluation of q-factor related quantities.
//! + [`Current`]: Evaluation of plasma current related quantities.
//! + [`Bfield`]: Evaluation of magnetic field related quantities.
//! + [`Mode`]: Single perturbation mode related quantities computation.
//!
//! ## Caching
//!
//! The trait [`Mode`]'s methods requires a [`ModeCache`] object to be passed as a parameter. Such an
//! object caches values such as angles' modulos and their sines/cosines, or amplitudes/phases
//! calculated with interpolation. It may also store the necessary
//! [`Accelerators`](rsl_interpolation::Accelerator). Since many evaluation methods are called with
//! the same arguments sequentially (like in the case of the perturbed equations of motion), it
//! makes sense to cache values that appear many times in these methods and can be expensive in
//! tight loops.
//!
//! + [`FluteModeCache`]: Cache for [`FluteMode`]
//! + [`NcFluteModeCache`]: Cache for [`NcFluteMode`]
//!
//! ## Data extraction
//!
//! The [`extract`] module provides methods for extracting data arrays and
//! [`Variables`](netcdf::Variable) from the netCDF file.
//!
//! + [`extract::open`]: Open a netCDF file.
//! + [`extract::scalar`]: Scalar value extraction.
//! + [`extract::array_1d`]: 1D array extraction.
//! + [`extract::array_2d`]: 2D array extraction.
//! + [`extract::array_3d`]: 3D array extraction.
//! + [`extract::mode_arrays`]: Extraction of the α and φ arrays of the {m,n} mode. Modes are indexed
//! by their mode numbers, rather than the logical index they appear on the data
//! arrays.
//! + [`extract::variable`]: Extraction of a variable as a [`Variable`](netcdf::Variable).
//! + [`extract::attribute`]: Extraction of a file's attribute as a String.
//! + [`extract::version`]: Extraction of a files convention [`Semantic Version`](https://semver.org/)

mod error;
mod eval;
mod objects;

// ============== Public API

pub mod extract;

pub use error::{EqError, EvalError, NcError};

pub use objects::{EquilibriumType, LastClosedFluxSurface};

pub use objects::nc_flux::FluxCoordinateState;

pub use eval::ModeCache;
pub use eval::{Bfield, Current, FluxCommute, Geometry, Mode, Qfactor};
pub use eval::{DynMode, DynModeCache};

pub use objects::geometries::LarGeometry;
pub use objects::geometries::NcGeometry;
pub use objects::geometries::NcGeometryBuilder;

pub use objects::qfactors::NcQfactor;
pub use objects::qfactors::NcQfactorBuilder;
pub use objects::qfactors::ParabolicQfactor;
pub use objects::qfactors::UnityQfactor;

pub use objects::currents::LarCurrent;
pub use objects::currents::NcCurrent;
pub use objects::currents::NcCurrentBuilder;

pub use objects::bfield::LarBfield;
pub use objects::bfield::NcBfield;
pub use objects::bfield::NcBfieldBuilder;

pub use objects::flute_mode::{FluteMode, FluteModeCache};
pub use objects::nc_flute_mode::{NcFluteMode, NcFluteModeBuilder, NcFluteModeCache, PhaseMethod};

pub use objects::perturbation::{DynModeCaches, DynModes, Perturbation};

// ============== Configuration constants

/// Crate configuration constants.
pub mod constants {
    /// The default `B` array `θ` padding width.
    pub const DEFAULT_THETA_PADDING_WIDTH: usize = 15;

    /// The index of the flux array, under which to switch to the analytical formula for the
    /// [`crate::NcFluteMode`] in order to ensure the correct `~sqrt(ψ)` behaviour near the axis.
    pub const NC_ANALYTICAL_THRESHOLD_INDEX: usize = 3;
}
