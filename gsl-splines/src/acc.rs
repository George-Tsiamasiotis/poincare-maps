use crate::RgslInterpAccel;

use std::cell::RefCell;
use std::rc::Rc;

/// Thin wrapper around `rgsl::InterpAccel`. This object can be mutable and shared across many
/// splines that evaluate at the same x point simultaneously.
pub struct Accelerator {
    pub(crate) gsl_iterp_accel: RgslInterpAccel,
}

impl Accelerator {
    /// Creates a new `Accelerator`.
    pub fn new() -> Rc<RefCell<Self>> {
        //Calls `gsl_interp_accel_alloc()`. The rest is taken care by `rgsl`.

        Rc::new(RefCell::new(Accelerator {
            gsl_iterp_accel: RgslInterpAccel::new(),
        }))
    }

    /// Resets the `Accelerator`'s cache and stats.
    pub(crate) fn reset(&mut self) {
        self.gsl_iterp_accel.reset();
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
