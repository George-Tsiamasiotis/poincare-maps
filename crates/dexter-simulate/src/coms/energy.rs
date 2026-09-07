//! Calculation of the Energy in a 2D grid.

#![expect(clippy::min_ident_chars, reason = "Hamiltonian terms")]

use ndarray::{Array1, Array2};
use rsl_interpolation::Accelerator2d;
use std::f64::consts::TAU;

use dexter_machine::Machine;

use crate::COMError;
use crate::MagneticFlux::*;
use crate::coms::COMs;

/// Calculation the Energy on a 2D meshgrid of the `ψ` and `θ` arrays, in Normalized Units.
///
/// Both `pzeta` and `mu` fields must be defined.
///
/// # Errors
///
/// Returns a [`COMError`] if `pzeta` or `mu` are not defined, or if any of the
/// `psi_array` values is out of bounds.
pub(crate) fn energy_of_psi_grid(
    coms: &COMs,
    machine: Machine,
    psi_array: &Array1<f64>,
    theta_array: &Array1<f64>,
) -> Result<Array2<f64>, COMError> {
    let Some(pzeta) = coms.pzeta else {
        return Err(COMError::UndefinedPzeta);
    };
    let Some(mu) = coms.mu else {
        return Err(COMError::UndefinedMu);
    };

    let mut grid = Array2::from_elem((psi_array.len(), theta_array.len()), f64::NAN);

    let acc = &mut Accelerator2d::new();

    let mod_theta_array = theta_array.mapv(|theta| theta.rem_euclid(TAU));

    // Iterate though `psi_array` first to avoid unnecessarily recalculating `rho`.
    // Use `<>_of_psi` evaluation methods to avoid error propagation through double interpolations.
    for i in 0..psi_array.len() {
        let psi = Toroidal(psi_array[i]);
        let psip = machine.qfactor().eval_other(psi, acc.xacc())?;
        let g = machine.current().eval_g(psi, acc.xacc())?;
        let rho = (pzeta + psip.value()) / g;
        for j in 0..mod_theta_array.len() {
            let theta = mod_theta_array[j];
            let b = machine.bfield().eval_b(psi, theta, acc)?;
            grid[[i, j]] = (rho * b).powi(2) / 2.0 + mu * b
        }
    }

    Ok(grid)
}

/// Calculates the Energy on a 2D meshgrid of the `ψp` and `θ` arrays, in Normalized Units.
///
/// Both `pzeta` and `mu` fields must be defined.
///
/// Note that the [`Qfactor`] object is not needed here.
///
/// # Errors
///
/// Returns a [`COMError`] if `pzeta` or `mu` are not defined, or if any of the
/// `psi_array` values is out of bounds.
pub(crate) fn energy_of_psip_grid(
    coms: &COMs,
    machine: Machine,
    psip_array: &Array1<f64>,
    theta_array: &Array1<f64>,
) -> Result<Array2<f64>, COMError> {
    let Some(pzeta) = coms.pzeta else {
        return Err(COMError::UndefinedPzeta);
    };
    let Some(mu) = coms.mu else {
        return Err(COMError::UndefinedMu);
    };

    let mut grid = Array2::from_elem((psip_array.len(), theta_array.len()), f64::NAN);

    let acc = &mut Accelerator2d::new();

    let mod_theta_array = theta_array.mapv(|theta| theta.rem_euclid(TAU));

    // Iterate though `psi_array` first to avoid unnecessarily recalculating `rho`.
    for i in 0..psip_array.len() {
        let psip = Poloidal(psip_array[i]);
        let g = machine.current().eval_g(psip, acc.xacc())?;
        let rho = (pzeta + psip.value()) / g;
        for j in 0..mod_theta_array.len() {
            let theta = mod_theta_array[j];
            let b = machine.bfield().eval_b(psip, theta, acc)?;
            grid[[i, j]] = (rho * b).powi(2) / 2.0 + mu * b
        }
    }

    Ok(grid)
}

#[cfg(test)]
mod test {
    use super::*;
    use dexter_machine::extract::TEST_NETCDF_PATH;
    use dexter_machine::*;
    use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};
    use ndarray::{arr1, arr2};

    #[test]
    fn gcmotion_check() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let psi_array = arr1(&vec![0.01, 0.02]);
        let theta_array = arr1(&vec![-1.0, 1.0]);

        let coms = COMs {
            energy: None,
            pzeta: Some(-0.03),
            mu: Some(1e-4),
        };

        let grid_of_psi = coms
            .energy_of_psi_grid(machine, &psi_array, &theta_array)
            .unwrap();

        let expected = arr2(&[
            [0.00026296256388989685, 0.00026296256388989685],
            [0.00012897176092872728, 0.00012897176092872728],
        ])
        .to_owned();

        assert!(grid_of_psi.relative_eq(&expected, 1e-15, 1e-20));
    }

    #[test]
    fn toroidal_poloidal_equivalence() {
        let path = std::path::PathBuf::from(TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Steffen).build().unwrap();
        let current = NcCurrentBuilder::new(&path, Steffen).build().unwrap();
        let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let acc = &mut Accelerator::new();
        let psi_array = arr1(&vec![0.01, 0.02]);
        let theta_array = arr1(&vec![-1.0, 1.0]);
        let psip_array = psi_array.mapv(|psi| {
            machine
                .qfactor()
                .eval_other(Toroidal(psi), acc)
                .unwrap()
                .value()
        });

        let coms = COMs {
            energy: None,
            pzeta: Some(-0.03),
            mu: Some(1e-4),
        };

        let grid_of_psi = coms
            .energy_of_psi_grid(machine, &psi_array, &theta_array)
            .unwrap();
        let grid_of_psip = coms
            .energy_of_psip_grid(machine, &psip_array, &theta_array)
            .unwrap();

        assert!(grid_of_psi.relative_eq(&grid_of_psip, 1e-11, 1e-14));
    }
}
