use config::PBAR_STYLE;
use equilibrium::{Bfield, Current, Perturbation, Qfactor};

use indicatif::{ParallelProgressIterator, ProgressStyle};
use particle::{Mapping, Particle};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::{PoincareInit, Result, results::PoincareResults};

#[non_exhaustive]
pub struct Poincare {
    pub init: PoincareInit,
    pub particles: Vec<Particle>,
    pub results: PoincareResults,
}

impl Poincare {
    pub fn new(init: PoincareInit) -> Self {
        let particles = init.to_particles();
        Self {
            init,
            particles,
            results: PoincareResults::default(),
        }
    }

    pub fn run(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
        mapping: &Mapping,
    ) -> Result<()> {
        self.particles
            .par_iter_mut()
            .progress_with_style(
                ProgressStyle::with_template(PBAR_STYLE).unwrap_or(ProgressStyle::default_bar()),
            )
            .try_for_each(|p| p.map(qfactor, bfield, current, per, mapping))?;

        self.results = PoincareResults::new(&self, mapping);
        Ok(())
    }
}
