use std::time::Duration;

use config::PBAR_STYLE;
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use utils::array2D_getter_impl;

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use particle::{MappingParameters, Particle};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::Result;
use crate::{Flux, Radians};
use crate::{PoincareInit, PoincareResults};

#[non_exhaustive]
/// Poincare map calculation.
pub struct Poincare {
    /// Initial conditions arrays.
    pub init: PoincareInit,
    /// Poincare map parameters.
    pub mapping: MappingParameters,
    /// Tracked [`Particle`]s.
    pub particles: Vec<Particle>,
    /// Integration results
    pub results: PoincareResults,
}

impl Poincare {
    /// Creates a new Poincare object, initializing all particles from the initial condtions
    /// arrays.
    pub fn new(init: PoincareInit, mapping: MappingParameters) -> Self {
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
        currents: &Currents,
        perturbation: &Perturbation,
    ) -> Result<()> {
        // `.progress_with()` seems to update the pbar **before** the map() method is called, so we
        // must create it and update it manually.
        let pbar = ProgressBar::new(self.particles.len() as u64).with_style(
            ProgressStyle::with_template(PBAR_STYLE).unwrap_or(ProgressStyle::default_bar()),
        );
        pbar.enable_steady_tick(Duration::from_millis(100));

        self.particles.par_iter_mut().try_for_each(|p| {
            p.map(qfactor, bfield, currents, perturbation, &self.mapping)
                .inspect(|()| pbar.inc(1))
        })?;

        self.results = PoincareResults::new(self, &self.mapping);
        Ok(())
    }
}

impl Poincare {
    array2D_getter_impl!(angles, results.angles, Radians);
    array2D_getter_impl!(fluxes, results.fluxes, Flux);
}
