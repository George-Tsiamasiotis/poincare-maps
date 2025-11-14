use std::time::Duration;

use ndarray::{Array1, Array2, Axis};
use particle::{MappingParameters, Particle, PoincareSection};
use safe_unwrap::safe_unwrap;
use utils::array2D_getter_impl;

use crate::Poincare;
use crate::Result;
use crate::{Flux, Radians};

/// Stores the results of the Poincare map calculation
#[derive(Clone)]
pub struct PoincareResults {
    /// The calculated θ angles, ignoring [`PoincareSection`].
    pub thetas: Array2<Radians>,
    /// The calculated ζ angles, ignoring [`PoincareSection`].
    pub zetas: Array2<Radians>,
    /// The calculated ψp flux, ignoring [`PoincareSection`].
    pub psips: Array2<Flux>,
    /// The calculated ψ flux, ignoring [`PoincareSection`].
    pub psis: Array2<Flux>,
    /// Poincare map parameters.
    pub params: MappingParameters,
    total: usize,
    mapped: usize,
    escaped: usize,
    timed_out: usize,
    invalid_intersections: usize,
    failed: usize,
    /// Duration of the slowest particle.
    slowest: MapDuration,
    /// Duration of the fastest particle.
    fastest: MapDuration,
}

impl PoincareResults {
    /// Creates a [`PoincareResults`] from an already calculated Poincare map.
    pub fn new(poincare: &Poincare, params: &MappingParameters) -> Result<Self> {
        let mut results = Self {
            params: *params,
            ..Default::default()
        };
        results.store_arrays(poincare)?;
        results.calculate_particle_nums(poincare);
        results.calculate_durations(poincare);
        Ok(results)
    }

    /// Extracts angles and fluxes data from a [`Poincare`] and stores then in `self`.
    pub fn store_arrays(&mut self, poincare: &Poincare) -> Result<()> {
        // We dont now how many particle's got completely integrated, so we push a new row for
        // every successful one.
        // We also include the initial point for now and drop it later, otherwise the code gets
        // ugly.
        let columns = self.params.intersections + 1;
        let shape = (0, columns);
        self.zetas = Array2::from_elem(shape, Radians::NAN);
        self.psips = Array2::from_elem(shape, Flux::NAN);
        self.thetas = Array2::from_elem(shape, Radians::NAN);
        self.psis = Array2::from_elem(shape, Flux::NAN);

        /// Copies the array of the calculated evolution `source` data into a new 1D array with length
        /// `columns` and pushes it to the 2D array `array`. If `len(source) < columns`, which will
        /// happen with escaped or timed out particles, the rest of the array is filled with NaN. This
        /// allows us to plot those particle's as well, while keeping all the data in the same 2D
        /// array.
        macro_rules! copy_and_fill_with_nan_and_push_row {
            ($particle:ident, $results_array:ident, $source:ident) => {
                assert!($particle.evolution.steps_stored() <= columns);
                self.$results_array.push_row(
                    Array1::from_shape_fn(columns, |i| {
                        $particle
                            .evolution
                            .$source()
                            .get(i)
                            .copied()
                            .unwrap_or(f64::NAN)
                    })
                    .view(),
                )?;
                // Still includes the initial point
                assert_eq!(self.$results_array.ncols(), self.params.intersections + 1);
            };
        }
        for p in poincare.particles.iter() {
            if should_be_plotted(p) {
                copy_and_fill_with_nan_and_push_row!(p, zetas, zeta);
                copy_and_fill_with_nan_and_push_row!(p, psips, psip);
                copy_and_fill_with_nan_and_push_row!(p, thetas, theta);
                copy_and_fill_with_nan_and_push_row!(p, psis, psi);
            }
        }

        // Remove intial points
        self.zetas.remove_index(Axis(1), 0);
        self.psips.remove_index(Axis(1), 0);
        self.thetas.remove_index(Axis(1), 0);
        self.psis.remove_index(Axis(1), 0);

        Ok(())
    }

    /// Counts the occurences of each [`IntegrationStatus`]'s variants.
    fn calculate_particle_nums(&mut self, poincare: &Poincare) {
        macro_rules! count_variants {
            ($is_enum:ident) => {
                poincare
                    .particles
                    .iter()
                    .filter(|p| p.status.$is_enum())
                    .count()
            };
        }
        self.mapped = count_variants!(is_mapped);
        self.escaped = count_variants!(is_escaped);
        self.timed_out = count_variants!(is_timed_out);
        self.escaped = count_variants!(is_escaped);
        self.invalid_intersections = count_variants!(is_invalid_intersections);
        self.failed = count_variants!(is_failed);
        self.total = poincare.particles.len();
    }

    fn calculate_durations(&mut self, poincare: &Poincare) {
        let slowest = safe_unwrap!(
            "poincare.particles is non-empty",
            poincare
                .particles
                .iter()
                .max_by_key(|p| p.evolution.duration)
        );
        let fastest = safe_unwrap!(
            "poincare.particles is non-empty",
            poincare
                .particles
                .iter()
                .filter(|p| p.evolution.steps_stored() > 0) // Drop invalid
                .min_by_key(|p| p.evolution.duration)
        );
        self.slowest = MapDuration::from(slowest);
        self.fastest = MapDuration::from(fastest);
    }
}

// Data extraction
impl PoincareResults {
    // Make them availiable to [`Poincare`]
    array2D_getter_impl!(zetas, zetas, Radians);
    array2D_getter_impl!(psips, psips, Radians);
    array2D_getter_impl!(thetas, thetas, Flux);
    array2D_getter_impl!(psis, psis, Flux);

    /// Returns rhe calculated map's angles, corresponding to either `ζ` or `θ`, depending on the
    /// [`PoincareSection`]
    pub fn angles(&self) -> Array2<Radians> {
        match self.params.section {
            PoincareSection::ConstTheta => self.zetas.clone(),
            PoincareSection::ConstZeta => self.thetas.clone(),
        }
    }

    /// Returns rhe calculated map's fluxes, corresponding to either `ψp` or `ψ`, depending on the
    /// [`PoincareSection`]
    pub fn fluxes(&self) -> Array2<Flux> {
        match self.params.section {
            PoincareSection::ConstTheta => self.psips.clone(),
            PoincareSection::ConstZeta => self.psis.clone(),
        }
    }
}

/// Returns true if the particle should be plotted in the final Poincare map.
fn should_be_plotted(particle: &Particle) -> bool {
    let status_ok =
        particle.status.is_mapped() | particle.status.is_timed_out() | particle.status.is_escaped();
    // Drop particles that were initialized outside the wall
    let length_ok = particle.evolution.steps_stored() > 1;

    status_ok && length_ok
}

/// Helper struct to display fastest and slowest particles
#[derive(Default, Clone)]
struct MapDuration {
    pub steps: usize,
    pub duration: Duration,
}

impl From<&Particle> for MapDuration {
    fn from(p: &Particle) -> Self {
        Self {
            steps: p.evolution.steps_taken(),
            duration: p.evolution.duration,
        }
    }
}

impl std::fmt::Debug for MapDuration {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "duration: {:?} ({} steps)", self.duration, self.steps)
    }
}

impl Default for PoincareResults {
    fn default() -> Self {
        Self {
            thetas: Default::default(),
            zetas: Default::default(),
            psips: Default::default(),
            psis: Default::default(),
            // Will be replaced
            params: MappingParameters::new(PoincareSection::ConstTheta, 0.0, 0),
            total: Default::default(),
            mapped: Default::default(),
            escaped: Default::default(),
            timed_out: Default::default(),
            invalid_intersections: Default::default(),
            failed: Default::default(),
            slowest: Default::default(),
            fastest: Default::default(),
        }
    }
}

impl std::fmt::Debug for PoincareResults {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PoincareResults")
            .field("total particles", &self.total)
            .field("integrated", &self.mapped)
            .field("escaped", &self.escaped)
            .field("timed_out", &self.timed_out)
            .field("invalid_intersections", &self.invalid_intersections)
            .field("failed", &self.failed)
            .field("slowest", &self.slowest)
            .field("fastest", &self.fastest)
            .finish()
    }
}
