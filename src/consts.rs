/// The initial capacity of the evolution time series Vecs.
pub const EVOLUTION_INIT_CAPACITY: usize = 2000;

/// The initial time step for the RKF45 adaptive step method. Should be small
/// enough to account for fast particles. The value is empirical.
pub const RKF45_FIRST_STEP: f64 = 1e-4;

// The maximum amount of steps a particle can make before terminating its integration.
pub const MAX_STEPS: usize = 1e9 as usize;
