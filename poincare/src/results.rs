use ndarray::{Array2, Axis};
use particle::{IntegrationStatus, MappingParameters};
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
    // TODO: statistics
}

impl PoincareResults {
    /// Creates a [`PoincareResults`] from an already calculated Poincare map.
    pub fn new(poincare: &Poincare, params: &MappingParameters) -> Result<Self> {
        // We dont now how many particle's got completely integrated, so we push a new row for
        // every successful one.
        // We also include the initial point for now and drop it later, otherwise the code gets
        // ugly.
        let shape = (0, params.intersections + 1);
        let mut angles: Array2<Radians> = Array2::from_elem(shape, Radians::NAN);
        let mut fluxes: Array2<Flux> = Array2::from_elem(shape, Flux::NAN);

        for p in poincare.particles.iter() {
            use particle::PoincareSection::*;
            if matches!(p.status, IntegrationStatus::Integrated) {
                match params.section {
                    ConstTheta => {
                        angles.push_row(p.evolution.zeta().view())?;
                        fluxes.push_row(p.evolution.psip().view())?;
                    }
                    ConstZeta => {
                        angles.push_row(p.evolution.theta().view())?;
                        fluxes.push_row(p.evolution.psi().view())?;
                    }
                }
            }
        }
        // Remove intial points
        angles.remove_index(Axis(1), 0);
        fluxes.remove_index(Axis(1), 0);

        assert!(!angles.is_any_nan(), "Poincare calculation returned NaN");
        assert!(!fluxes.is_any_nan(), "Poincare calculation returned NaN");

        Ok(Self { angles, fluxes })
    }
}

impl PoincareResults {
    // Make them availiable to [`Poincare`]
    array2D_getter_impl!(angles, angles, Radians);
    array2D_getter_impl!(fluxes, fluxes, Flux);
}
