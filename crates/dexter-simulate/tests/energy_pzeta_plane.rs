//! Test the creation of an `EnergyPzetaPlane` from a set of `COMs`.

#![allow(unused_variables)]

use dexter_equilibrium::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
use dexter_equilibrium::*;
use dexter_simulate::*;
use parabola::Parabola;
use std::path::PathBuf;

#[test]
fn lar_energy_pzeta_parabola() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(UnityQfactor::new(lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::zero(),
    };

    let coms = COMs {
        energy: None,
        pzeta: Some(-0.02),
        mu: Some(7e-6),
    };

    let plane = coms.build_energy_pzeta_plane(&equilibrium).unwrap();
    let axis: &Parabola = plane.axis_parabola();
    let left_wall: &Parabola = plane.left_wall_parabola();
    let right_wall: &Parabola = plane.right_wall_parabola();
    let tp_boundary: &TrappedPassingBoundary = plane.tp_boundary();

    assert_eq!(plane.mu(), coms.mu.unwrap());
}

#[test]
fn toroidal_nc_energy_pzeta_parabola() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, "steffen").build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, "steffen").build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, "bicubic").build().unwrap()),
        perturbation: Perturbation::zero(),
    };

    let coms = COMs {
        energy: None,
        pzeta: Some(-0.02),
        mu: Some(7e-6),
    };

    let plane = coms.build_energy_pzeta_plane(&equilibrium).unwrap();
    let axis: &Parabola = plane.axis_parabola();
    let left_wall: &Parabola = plane.left_wall_parabola();
    let right_wall: &Parabola = plane.right_wall_parabola();
    let tp_boundary: &TrappedPassingBoundary = plane.tp_boundary();

    assert_eq!(plane.mu(), coms.mu.unwrap());
}

#[test]
fn poloidal_nc_energy_pzeta_parabola() {
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, "steffen").build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, "steffen").build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, "bicubic").build().unwrap()),
        perturbation: Perturbation::zero(),
    };

    let coms = COMs {
        energy: None,
        pzeta: Some(-0.02),
        mu: Some(7e-6),
    };

    let plane = coms.build_energy_pzeta_plane(&equilibrium).unwrap();
    let axis: &Parabola = plane.axis_parabola();
    let left_wall: &Parabola = plane.left_wall_parabola();
    let right_wall: &Parabola = plane.right_wall_parabola();
    let tp_boundary: &TrappedPassingBoundary = plane.tp_boundary();

    assert_eq!(plane.mu(), coms.mu.unwrap());
}
