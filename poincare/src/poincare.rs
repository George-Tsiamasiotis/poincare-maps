use std::time::Duration;

use config::PBAR_STYLE;
use equilibrium::{Bfield, Current, Perturbation, Qfactor};
use utils::array2D_getter_impl;

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use particle::{Mapping, Particle};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::{PoincareInit, Result, results::PoincareResults};

#[non_exhaustive]
/// Poincare map calculation.
pub struct Poincare {
    /// Initial conditions arrays.
    pub init: PoincareInit,
    /// Poincare map parameters.
    pub mapping: Mapping,
    /// Tracked [`Particle`]s.
    pub particles: Vec<Particle>,
    /// Integration results
    pub results: PoincareResults,
}

impl Poincare {
    /// Creates a new Poincare object, initializing all particles from the initial condtions
    /// arrays.
    pub fn new(init: PoincareInit, mapping: Mapping) -> Self {
        let particles = init.to_particles();
        Self {
            init,
            mapping,
            particles,
            results: PoincareResults::default(),
        }
    }

    /// Performs the calculation
    pub fn run(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> Result<()> {
        // `.progress_with()` seems to update the pbar **before** the map() method is called, so we
        // must create it and update it manually.
        let pbar = ProgressBar::new(self.particles.len() as u64).with_style(
            ProgressStyle::with_template(PBAR_STYLE).unwrap_or(ProgressStyle::default_bar()),
        );
        pbar.enable_steady_tick(Duration::from_millis(100));

        self.particles.par_iter_mut().try_for_each(|p| {
            p.map(qfactor, bfield, current, per, &self.mapping)
                .inspect(|()| pbar.inc(1))
        })?;

        self.results = PoincareResults::new(self, &self.mapping);
        Ok(())
    }

    pub fn results(&self) -> PoincareResults {
        self.results.clone()
    }
}

array2D_getter_impl!(Poincare, angles, results.angles);
array2D_getter_impl!(Poincare, fluxes, results.fluxes);
