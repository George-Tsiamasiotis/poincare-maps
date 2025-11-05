// ==================== Solver

/// The maximum amount of steps a particle can make before terminating its integration.
pub const MAX_STEPS: usize = 10_000_000;

/// The initial time step for the RKF45 adaptive step method. Should be small enough to account
/// for fast particles. The value is empirical.
pub const RKF45_FIRST_STEP: f64 = 1e-3;

/// The safety factor of the solver. Should be less than 1.0
pub const SAFETY_FACTOR: f64 = 0.9;

/// The relative tolerance of the energy error in every step.
pub const ENERGY_REL_TOL: f64 = 1e-11;

/// The relative tolerance of the stepping error in every step.
pub const STEP_REL_TOL: f64 = 1e-11;

// ==================== Mapping

/// The maximum allowed absolute difference between two consecutive intersections.
pub const MAP_THRESHOLD: f64 = 1e-9;

// ==================== Misc

/// The starting capacity of the Evolution time series vectors.
pub const EVOLUTION_INIT_CAPACITY: usize = 2000;

// ==================== Display

/// The Poincare map calculation progress bar style.
pub const PBAR_STYLE: &str =
    "[{elapsed_precise}] {wide_bar:.cyan/blue} {spinner} {pos:>4}/{len:4} {msg}";
