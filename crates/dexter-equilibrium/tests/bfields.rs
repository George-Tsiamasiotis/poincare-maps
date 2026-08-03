//! Test Bfields functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_equilibrium::*;
use ndarray::{Array1, Array2};

#[test]
#[rustfmt::skip]
fn lar_bfield() {
    let bfield = dbg!(LarBfield::new());
    assert_eq!(bfield.psi_state(), FluxCoordinateState::Good);
    assert_eq!(bfield.psip_state(), FluxCoordinateState::Bad);

    let acc = &mut Accelerator2d::new();
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = bfield.b_of_psi(psi, theta, acc).unwrap();
    let _: f64 = bfield.db_dpsi(psi, theta, acc).unwrap();
    let _: f64 = bfield.db_of_psi_dtheta(psi, theta, acc).unwrap();

    assert!(matches!(
        bfield.b_of_psip(psip, theta, acc),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        bfield.db_dpsip(psip, theta, acc),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        bfield.db_of_psip_dtheta(psip, theta, acc),
        Err(EvalError::UndefinedEvaluation(..))
    ));
}

#[test]
fn nc_bfield_no_pad() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let interp_type = Interpolation2dType::Bicubic;
    let builder = NcBfieldBuilder::new(&path, interp_type).with_padding(0);
    let bfield = dbg!(builder.build().unwrap());
    assert_eq!(bfield.psi_state(), FluxCoordinateState::Good);
    assert_eq!(bfield.psip_state(), FluxCoordinateState::Good);

    let equilibrium_type: ObjectType = bfield.object_type();
    let netcdf_version: semver::Version = bfield.netcdf_version();
    let path: PathBuf = bfield.path();
    let interp_type: Interpolation2dType = bfield.interp_type();
    let baxis: f64 = bfield.baxis();
    let shape: (usize, usize) = bfield.shape();
    let psi_state: FluxCoordinateState = bfield.psi_state();
    let psip_state: FluxCoordinateState = bfield.psip_state();
    let psi_last: f64 = bfield.psi_last().unwrap();
    let psip_last: f64 = bfield.psip_last().unwrap();
    let psi_array: Array1<f64> = bfield.psi_array().unwrap();
    let psip_array: Array1<f64> = bfield.psip_array().unwrap();
    let theta_array: Array1<f64> = bfield.theta_array();
    let theta_array_padded: Array1<f64> = bfield.theta_array_padded();
    let b_array: Array2<f64> = bfield.b_array();
    let b_array_padded: Array2<f64> = bfield.b_array_padded();

    assert_eq!(theta_array, theta_array_padded);
    assert_eq!(b_array, b_array_padded);

    let acc = &mut Accelerator2d::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = bfield.b_of_psi(psi, theta, acc).unwrap();
    let _: f64 = bfield.b_of_psip(psip, theta, acc).unwrap();
    let _: f64 = bfield.db_dpsi(psi, theta, acc).unwrap();
    let _: f64 = bfield.db_dpsip(psip, theta, acc).unwrap();
    let _: f64 = bfield.db_of_psi_dtheta(psi, theta, acc).unwrap();
    let _: f64 = bfield.db_of_psip_dtheta(psip, theta, acc).unwrap();
}

#[test]
fn nc_bfield_pad() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let interp_type = Interpolation2dType::Bicubic;
    let no_pad_builder = NcBfieldBuilder::new(&path, interp_type).with_padding(0);
    let pad_builder = NcBfieldBuilder::new(&path, interp_type).with_padding(10);
    let no_pad_bfield = dbg!(no_pad_builder.build().unwrap());
    let pad_bfield = dbg!(pad_builder.build().unwrap());

    assert_eq!(no_pad_bfield.padding(), 0);
    assert_eq!(pad_bfield.padding(), 10);

    let netcdf_shape = no_pad_bfield.shape();

    assert_eq!(netcdf_shape.0, pad_bfield.shape().0);
    assert_eq!(
        netcdf_shape.1,
        pad_bfield.shape_padded().1 - 2 * pad_bfield.padding()
    );
}
