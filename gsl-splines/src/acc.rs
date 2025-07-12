use crate::RgslInterpAccel;

/// Thin wrapper around `rgsl::InterpAccel`. This object can be mutable and shared across many
/// splines that evaluate at the same x point simultaneously.
pub struct Accelerator {
    pub(crate) gsl_iterp_accel: RgslInterpAccel,
}

impl Accelerator {
    /// Creates a new `Accelerator`
    pub fn new() -> Self {
        //Calls `gsl_interp_accel_alloc()`. The rest is taken care by `rgsl`.
        Accelerator {
            gsl_iterp_accel: RgslInterpAccel::new(),
        }
    }

    /// Resets the `Accelerator`'s cache and stats.
    pub fn reset(&mut self) {
        self.gsl_iterp_accel.reset();
    }
}

impl Default for Accelerator {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for Accelerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Accelerator")
            .field("hits  ", &self.gsl_iterp_accel.0.hit_count)
            .field("misses", &self.gsl_iterp_accel.0.miss_count)
            .finish()
    }
}
