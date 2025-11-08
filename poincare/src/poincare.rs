use std::time::Duration;

use config::{POINCARE_PBAR_STYLE, POINCARE_PROGRESS_CHARS};
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use particle::{MappingParameters, Particle};
use utils::array2D_getter_impl;

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
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
    pub params: MappingParameters,
    /// Tracked [`Particle`]s.
    pub particles: Vec<Particle>,
    /// Integration results
    pub results: PoincareResults,
}

impl Poincare {
    /// Creates a new Poincare object, initializing all particles from the initial condtions
    /// arrays.
    pub fn new(init: PoincareInit, params: MappingParameters) -> Self {
        let particles = init.to_particles();
        Self {
            init,
            params,
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
            ProgressStyle::with_template(POINCARE_PBAR_STYLE)
                .unwrap_or(ProgressStyle::default_bar())
                .progress_chars(POINCARE_PROGRESS_CHARS),
        );
        pbar.enable_steady_tick(Duration::from_millis(100));
        pbar.println(format!(
            "🗿 Integrating {} particles for {} intersections",
            self.particles.len(),
            self.params.intersections
        ));
        pbar.println(format!("🚀 Using {} threads", rayon::current_num_threads()));

        self.particles.par_iter_mut().try_for_each(|p| {
            p.map(qfactor, bfield, currents, perturbation, &self.params)
                .inspect(|()| pbar.inc(1))
        })?;
        pbar.println("✅️ Integration Done");
        pbar.finish();

        self.results = PoincareResults::new(self, &self.params)?;
        Ok(())
    }
}

impl Poincare {
    array2D_getter_impl!(angles, results.angles, Radians);
    array2D_getter_impl!(fluxes, results.fluxes, Flux);
}

impl std::fmt::Debug for Poincare {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Poincare")
            .field("params", &self.params)
            .field("results", &self.results)
            .finish()
    }
}
