//! Representation of an equilibrium and added perturbations.

use crate::{Bfield, Current, Geometry, Perturbation, Qfactor};

/// An Equilibrium, containing all information about the geometry, currents, qfactor and magnetic
/// field of the configuration, as well as any present perturbation.
#[derive(Debug)]
pub struct Equilibrium {
    /// The Equilibrium's [`Geometry`].
    pub geometry: Option<Box<dyn Geometry>>,
    /// The Equilibrium's [`Qfactor`].
    pub qfactor: Box<dyn Qfactor>,
    /// The Equilibrium's [`Current`].
    pub current: Box<dyn Current>,
    /// The Equilibrium's [`Bfield`].
    pub bfield: Box<dyn Bfield>,
    /// The Equilibrium's [`Perturbation`].
    pub perturbation: Perturbation,
}

impl Equilibrium {
    /// Returns the value of the last closed toroidal flux `ψ_last`.
    #[must_use]
    pub fn psi_last(&self) -> f64 {
        self.qfactor.psi_last()
    }

    /// Returns the value of the last closed toroidal flux `ψp_last`.
    #[must_use]
    pub fn psip_last(&self) -> f64 {
        self.qfactor.psip_last()
    }
}
