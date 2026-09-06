//! Test Currents functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use approx::assert_abs_diff_eq;
use dexter_machine::*;
use ndarray::Array1;

#[test]
fn lar_current() {
    let current = dbg!(LarCurrent::new());
    assert_eq!(current.psi_state(), FluxCoordinateState::Good);
    assert_eq!(current.psip_state(), FluxCoordinateState::Good);

    let mut acc = Accelerator::new();
    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert_abs_diff_eq!(current.eval_g(psi, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(current.eval_g(psip, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(current.eval_i(psi, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_g_deriv(psi, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_g_deriv(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i_deriv(psi, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i_deriv(psip, &mut acc).unwrap(), 0.0);
}

#[test]
fn nc_current() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Cubic;
    let builder = NcCurrentBuilder::new(&path, typ);
    let current = dbg!(builder.build().unwrap());
    assert_eq!(current.psi_state(), FluxCoordinateState::Good);
    assert_eq!(current.psip_state(), FluxCoordinateState::Good);

    let machine_type: MachineType = current.machine_type();
    let netcdf_version: semver::Version = current.netcdf_version();
    let path: PathBuf = current.path();
    let interp_type: Interpolation1dType = current.interp_type();
    let psi_state: FluxCoordinateState = current.psi_state();
    let psip_state: FluxCoordinateState = current.psip_state();
    let psi_last: f64 = current.psi_last().unwrap();
    let psip_last: f64 = current.psip_last().unwrap();
    let psi_array: Array1<f64> = current.psi_array().unwrap();
    let psip_array: Array1<f64> = current.psip_array().unwrap();
    let g_array: Array1<f64> = current.g_array();
    let i_array: Array1<f64> = current.i_array();

    let mut acc = Accelerator::new();
    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    let _: f64 = current.eval_g(psi, &mut acc).unwrap();
    let _: f64 = current.eval_g(psip, &mut acc).unwrap();
    let _: f64 = current.eval_i(psi, &mut acc).unwrap();
    let _: f64 = current.eval_i(psip, &mut acc).unwrap();
    let _: f64 = current.eval_i_deriv(psi, &mut acc).unwrap();
    let _: f64 = current.eval_i_deriv(psip, &mut acc).unwrap();
    let _: f64 = current.eval_i_deriv(psi, &mut acc).unwrap();
    let _: f64 = current.eval_i_deriv(psip, &mut acc).unwrap();
}
