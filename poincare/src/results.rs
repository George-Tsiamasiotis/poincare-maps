use ndarray::Array1;
use ndarray::{Array2, Axis};
use particle::{MappingParameters, Particle, PoincareSection};
use utils::array2D_getter_impl;

use crate::Poincare;
use crate::Result;
use crate::{Flux, Radians};

/// Stores the results of the Poincare map calculation
#[derive(Default, Clone)]
pub struct PoincareResults {
    /// The calculated angles, corresponding to either `zeta` or `theta`, depending on the
    /// [`PoincareSection`]
    pub angles: Array2<Radians>,
    /// The calculated fluxes, corresponding to either `psip` or `psi`, depending on the
    /// [`PoincareSection`]
    pub fluxes: Array2<Flux>,
    total: usize,
    integrated: usize,
    escaped: usize,
    timed_out: usize,
    invalid_intersections: usize,
    failed: usize,
}

impl PoincareResults {
    /// Creates a [`PoincareResults`] from an already calculated Poincare map.
    pub fn new(poincare: &Poincare, params: &MappingParameters) -> Result<Self> {
        let (angles, fluxes) = calculate_arrays(poincare, params)?;
        let mut results = Self {
            angles,
            fluxes,
            ..Default::default()
        };
        results.calculate_particle_nums(poincare);
        Ok(results)
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
        self.integrated = count_variants!(is_integrated);
        self.escaped = count_variants!(is_escaped);
        self.timed_out = count_variants!(is_timed_out);
        self.escaped = count_variants!(is_escaped);
        self.invalid_intersections = count_variants!(is_invalid_intersections);
        self.failed = count_variants!(is_failed);
        self.total = poincare.particles.len();
    }
}

impl PoincareResults {
    // Make them availiable to [`Poincare`]
    array2D_getter_impl!(angles, angles, Radians);
    array2D_getter_impl!(fluxes, fluxes, Flux);
}

/// Extracts angle and flux data from a [`Poincare`] object and returns the 2 2D arrays.
pub fn calculate_arrays(
    poincare: &Poincare,
    params: &MappingParameters,
) -> Result<(Array2<f64>, Array2<f64>)> {
    // We dont now how many particle's got completely integrated, so we push a new row for
    // every successful one.
    // We also include the initial point for now and drop it later, otherwise the code gets
    // ugly.
    let columns = params.intersections + 1;
    let shape = (0, columns);
    let mut angles: Array2<Radians> = Array2::from_elem(shape, Radians::NAN);
    let mut fluxes: Array2<Flux> = Array2::from_elem(shape, Flux::NAN);

    /// Copies the array of the calculated evolution `source` data into a new 1D array with length
    /// `columns` and pushes it to the 2D array `array`. If `len(source) < columns`, which will
    /// happen with escaped or timed out particles, the rest of the array is filled with NaN. This
    /// allows us to plot those particle's as well, while keeping all the data in the same 2D
    /// array.
    macro_rules! copy_and_fill_with_nan_and_push_row {
        ($particle:ident, $results_array:ident, $source:ident) => {
            assert!($particle.evolution.steps_stored() <= columns);
            $results_array.push_row(
                Array1::from_shape_fn(columns, |i| {
                    $particle
                        .evolution
                        .$source()
                        .get(i)
                        .copied()
                        .unwrap_or(f64::NAN)
                        .clone()
                })
                .view(),
            )?;
        };
    }

    for p in poincare.particles.iter() {
        if should_be_plotted(p) {
            match params.section {
                PoincareSection::ConstTheta => {
                    copy_and_fill_with_nan_and_push_row!(p, angles, zeta);
                    copy_and_fill_with_nan_and_push_row!(p, fluxes, psip);
                }
                PoincareSection::ConstZeta => {
                    copy_and_fill_with_nan_and_push_row!(p, angles, theta);
                    copy_and_fill_with_nan_and_push_row!(p, fluxes, psi);
                }
            }
        }
    }
    // Remove intial points
    angles.remove_index(Axis(1), 0);
    fluxes.remove_index(Axis(1), 0);

    assert_eq!(angles.ncols(), params.intersections);
    assert_eq!(angles.shape(), fluxes.shape());

    Ok((angles, fluxes))
}

/// Returns true if the particle should be plotted in the final Poincare map.
fn should_be_plotted(particle: &Particle) -> bool {
    let status_ok = particle.status.is_integrated()
        | particle.status.is_timed_out()
        | particle.status.is_escaped();
    // Drop particles that were initialized outside the wall
    let length_ok = particle.evolution.steps_stored() > 1;

    status_ok && length_ok
}

impl std::fmt::Debug for PoincareResults {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PoincareResults")
            .field("total particles", &self.total)
            .field("integrated", &self.integrated)
            .field("escaped", &self.escaped)
            .field("timed_out", &self.timed_out)
            .field("invalid_intersections", &self.invalid_intersections)
            .field("failed", &self.failed)
            .finish()
    }
}
