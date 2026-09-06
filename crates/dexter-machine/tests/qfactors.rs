//! Test Qfactors functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use approx::{assert_abs_diff_eq, assert_relative_eq};
use dexter_machine::{MagneticFlux::*, *};
use ndarray::Array1;

#[test]
fn unity_qfactor() {
    let qfactor = dbg!(UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.45)));

    assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
    assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

    let mut acc = Accelerator::new();
    let psi = Toroidal(0.01);
    let psip = Poloidal(0.015);
    assert_abs_diff_eq!(qfactor.eval_q(psi, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_q(psip, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_other(psi, &mut acc).unwrap(), Poloidal(0.01));
    assert_abs_diff_eq!(qfactor.eval_other(psip, &mut acc).unwrap(), Toroidal(0.015));
    assert_abs_diff_eq!(qfactor.eval_deriv_of_other(psi, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_deriv_of_other(psip, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_iota(psi, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_iota(psip, &mut acc).unwrap(), 1.0);

    assert!(qfactor.eval_psi_of_q(1.0, &mut acc).is_err());
    assert!(qfactor.eval_psip_of_q(1.0, &mut acc).is_err());

    assert_eq!(qfactor.psi_last(), Toroidal(0.45));
    assert_eq!(qfactor.psip_last(), Poloidal(0.45));
}

#[test]
fn parabolic_qfactor() {
    let qfactor = dbg!(ParabolicQfactor::new(
        1.1,
        3.8,
        LastClosedFluxSurface::Toroidal(0.45)
    ));

    assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
    assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

    let qaxis: f64 = qfactor.qaxis();
    let qlast: f64 = qfactor.qlast();
    let psi_last: MagneticFlux = qfactor.psi_last();
    let psip_last: MagneticFlux = qfactor.psip_last();

    let mut acc = Accelerator::new();
    let psi = Toroidal(0.01);
    let psip = Poloidal(0.015);
    let _: f64 = qfactor.eval_q(psi, &mut acc).unwrap();
    let _: f64 = qfactor.eval_q(psip, &mut acc).unwrap();
    let _: f64 = qfactor.eval_iota(psi, &mut acc).unwrap();
    let _: f64 = qfactor.eval_iota(psip, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_other(psi, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_other(psip, &mut acc).unwrap();
    let _: f64 = qfactor.eval_deriv_of_other(psi, &mut acc).unwrap();
    let _: f64 = qfactor.eval_deriv_of_other(psip, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_psi_of_q(1.3, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_psip_of_q(1.3, &mut acc).unwrap();
}

#[test]
fn nc_qfactor() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Cubic;
    let builder = NcQfactorBuilder::new(&path, typ);
    let qfactor = dbg!(builder.build().unwrap());

    assert_eq!(qfactor.psi_state(), FluxCoordinateState::Good);
    assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

    let machine_type: MachineType = qfactor.machine_type();
    let netcdf_version: semver::Version = qfactor.netcdf_version();
    let path: PathBuf = qfactor.path();
    let interp_type: Interpolation1dType = qfactor.interp_type();
    let psi_state: FluxCoordinateState = qfactor.psi_state();
    let psip_state: FluxCoordinateState = qfactor.psip_state();
    let qaxis: f64 = qfactor.qaxis();
    let qlast: f64 = qfactor.qlast();
    let psi_last: MagneticFlux = qfactor.psi_last();
    let psip_last: MagneticFlux = qfactor.psip_last();
    let psi_array: Array1<f64> = qfactor.psi_array().unwrap();
    let psip_array: Array1<f64> = qfactor.psip_array().unwrap();
    let q_array: Array1<f64> = qfactor.q_array();

    let mut acc = Accelerator::new();
    let psi = Toroidal(0.01);
    let psip = Poloidal(0.015);
    let _: f64 = qfactor.eval_q(psi, &mut acc).unwrap();
    let _: f64 = qfactor.eval_q(psip, &mut acc).unwrap();
    let _: f64 = qfactor.eval_iota(psi, &mut acc).unwrap();
    let _: f64 = qfactor.eval_iota(psip, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_other(psi, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_other(psip, &mut acc).unwrap();
    let _: f64 = qfactor.eval_deriv_of_other(psi, &mut acc).unwrap();
    let _: f64 = qfactor.eval_deriv_of_other(psip, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_psi_of_q(1.3, &mut acc).unwrap();
    let _: MagneticFlux = qfactor.eval_psip_of_q(1.3, &mut acc).unwrap();
}

#[test]
fn nc_qfactor_inverse_q() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Akima;
    let builder = NcQfactorBuilder::new(&path, typ);
    let qfactor = dbg!(builder.build().unwrap());
    let acc = &mut Accelerator::new();

    let psi = Toroidal(0.1);
    let q_of_psi = qfactor.eval_q(psi, acc).unwrap();
    let psi_inverse = qfactor.eval_psi_of_q(q_of_psi, acc).unwrap();
    assert_relative_eq!(psi, psi_inverse, epsilon = 1e-6);

    let psip = Poloidal(0.1);
    let q_of_psip = qfactor.eval_q(psip, acc).unwrap();
    let psip_inverse = qfactor.eval_psip_of_q(q_of_psip, acc).unwrap();
    assert_relative_eq!(psip, psip_inverse, epsilon = 1e-6);
}
