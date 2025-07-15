use crate::RgslInterpAccel;

use std::cell::RefCell;
use std::rc::Rc;

/// An object that caches the previous interpolation lookup.
///
/// Splines that tend to evaluate around the same point frequently can get a significant
/// performance boost (up to x10, depending on the number of points) with the use of an
/// Accelerator.
///
/// This object is a thin wrapper around [`rgsl's InterpAccel`]. It is mutable and can be shared across
/// many splines defined over the same data points and evaluate at the same x point simultaneously.
///
/// ## Example
///
/// See [`example`]
///
/// [`example`]: ./index.html
/// [`rgsl's InterpAccel`]: https://docs.rs/GSL/latest/rgsl/types/interpolation/struct.InterpAccel.html

pub struct Accelerator {
    pub(crate) gsl_iterp_accel: RgslInterpAccel,
}

impl Accelerator {
    /// Creates a new Accelerator.
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_accelarator() {
        Accelerator::new();
    }
}
