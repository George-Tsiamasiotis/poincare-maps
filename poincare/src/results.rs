use ndarray::{Array2, Axis};
use particle::Mapping;
use utils::array2D_getter_impl;

use crate::Poincare;

/// Stores the results of the Poincare map calculation
#[derive(Default, Clone)]
pub struct PoincareResults {
    /// The calculated angles, corresponding to either `zeta` or `theta`, depending on the
    /// [`PoincareSection`]
    pub angles: Array2<f64>,
    /// The calculated fluxes, corresponding to either `psip` or `psi`, depending on the
    /// [`PoincareSection`]
    pub fluxes: Array2<f64>,
    // TODO: statistics
}

impl PoincareResults {
    /// Creates a [`PoincareResults`] from an already calculated Poincare map.
    pub fn new(poincare: &Poincare, mapping: &Mapping) -> Self {
        // Include initiali point
        let shape = (poincare.particles.len(), mapping.intersections + 1);
        let mut angles: Array2<f64> = Array2::zeros(shape);
        let mut fluxes: Array2<f64> = Array2::zeros(shape);

        poincare.particles.iter().enumerate().for_each(|(row, p)| {
            use particle::PoincareSection::*;
            match mapping.section {
                ConstTheta => {
                    angles.row_mut(row).assign(&p.evolution.zeta());
                    fluxes.row_mut(row).assign(&p.evolution.psip());
                }
                ConstZeta => {
                    angles.row_mut(row).assign(&p.evolution.theta());
                    fluxes.row_mut(row).assign(&p.evolution.psi());
                }
            }
        });
        // Remove intial points
        angles.remove_index(Axis(1), 0);
        fluxes.remove_index(Axis(1), 0);
        Self { angles, fluxes }
    }
}

// Make them availiable to Poincare
array2D_getter_impl!(PoincareResults, angles, angles);
array2D_getter_impl!(PoincareResults, fluxes, fluxes);
