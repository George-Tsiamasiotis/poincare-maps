use ndarray::{Array2, Axis};
use particle::Mapping;

use crate::Poincare;

#[derive(Default)]
pub struct PoincareResults {
    pub angles: Array2<f64>,
    pub fluxes: Array2<f64>,
}

impl PoincareResults {
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
        debug_assert_eq!(
            angles.shape(),
            &[poincare.particles.len(), mapping.intersections]
        );
        Self { angles, fluxes }
    }
}
