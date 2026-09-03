//! Calculations on the Constants of Motion (COMs) space.

#![expect(dead_code, reason = "to be re-written")]

mod energy;
mod energy_pzeta_plane;
mod tp_boundary;

use ndarray::{Array1, Array2};

use dexter_machine::Machine;

pub(crate) use energy_pzeta_plane::EnergyPzetaPlane;
pub(crate) use tp_boundary::TrappedPassingBoundary;

use crate::COMError;

/// The constants of motion `(E, Pζ, μ)` in an unperturbed equilibrium.
#[derive(Debug, Clone, Copy)]
pub(crate) struct COMs {
    /// The Energy in Normalized Units.
    pub energy: Option<f64>,
    /// The canonical momentum `Pζ` in Normalized Units.
    pub pzeta: Option<f64>,
    /// The magnetic moment `μ` in Normalized Units.
    pub mu: Option<f64>,
}

impl COMs {
    /// Calculates the Energy on a 2D meshgrid of the `ψ` and `θ` arrays, in Normalized Units.
    ///
    /// Both `pzeta` and `mu` fields must be defined.
    pub(crate) fn energy_of_psi_grid(
        &self,
        objects: Machine,
        psi_array: &Array1<f64>,
        theta_array: &Array1<f64>,
    ) -> Result<Array2<f64>, COMError> {
        energy::energy_of_psi_grid(self, objects, psi_array, theta_array)
    }

    /// Calculates the Energy on a 2D meshgrid of the `ψp` and `θ` arrays, in Normalized Units.
    ///
    /// Both `pzeta` and `mu` fields must be defined.
    ///
    /// Note that this calculation is independent from the [`Qfactor`](dexter_machine::Qfactor) used.
    pub(crate) fn energy_of_psip_grid(
        &self,
        objects: Machine,
        psip_array: &Array1<f64>,
        theta_array: &Array1<f64>,
    ) -> Result<Array2<f64>, COMError> {
        energy::energy_of_psip_grid(self, objects, psip_array, theta_array)
    }

    /// Constructs a [`EnergyPzetaPlane`].
    ///
    /// Only the [`COMs::mu`] field is necessary.
    pub(crate) fn build_energy_pzeta_plane(
        &self,
        objects: Machine,
    ) -> Result<EnergyPzetaPlane, COMError> {
        EnergyPzetaPlane::from_coms(objects, self)
    }
}
