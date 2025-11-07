use ndarray::{Array2, Axis};
use particle::{IntegrationStatus, MappingParameters};
use utils::array2D_getter_impl;

use crate::Poincare;
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
    // TODO: statistics
}

impl PoincareResults {
    /// Creates a [`PoincareResults`] from an already calculated Poincare map.
    pub fn new(poincare: &Poincare, params: &MappingParameters) -> Self {
        // Include initial point
        let shape = (poincare.particles.len(), params.intersections + 1);
        let mut angles: Array2<Radians> = Array2::from_elem(shape, Radians::NAN);
        let mut fluxes: Array2<Flux> = Array2::from_elem(shape, Flux::NAN);

        poincare.particles.iter().enumerate().for_each(|(row, p)| {
            use particle::PoincareSection::*;
            if matches!(p.status, IntegrationStatus::Integrated) {
                match params.section {
                    ConstTheta => {
                        angles.row_mut(row).assign(&p.evolution.zeta());
                        fluxes.row_mut(row).assign(&p.evolution.psip());
                    }
                    ConstZeta => {
                        angles.row_mut(row).assign(&p.evolution.theta());
                        fluxes.row_mut(row).assign(&p.evolution.psi());
                    }
                }
            }
        });
        // Remove intial points
        angles.remove_index(Axis(1), 0);
        fluxes.remove_index(Axis(1), 0);

        debug_assert!(!angles.is_any_nan(), "Poincare calculation returned NaN");
        debug_assert!(!fluxes.is_any_nan(), "Poincare calculation returned NaN");

        Self { angles, fluxes }
    }
}

impl PoincareResults {
    // Make them availiable to [`Poincare`]
    array2D_getter_impl!(angles, angles, Radians);
    array2D_getter_impl!(fluxes, fluxes, Flux);
}
