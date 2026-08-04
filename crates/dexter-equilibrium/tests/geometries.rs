//! Test Geometry functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_equilibrium::*;
use ndarray::{Array1, Array2};

#[test]
fn lar_geometry() {
    let geometry = dbg!(LarGeometry::new(2.0, 1.75, 0.5));

    assert_eq!(geometry.psi_state(), FluxCoordinateState::Good);
    assert_eq!(geometry.psip_state(), FluxCoordinateState::Bad);

    let equilibrium_type: ObjectType = geometry.object_type();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let rlast: f64 = geometry.rlast();
    let psi_last: f64 = geometry.psi_last().unwrap();
    assert!(geometry.psip_last().is_none());
    let rlab_last: Array1<f64> = geometry.rlab_last();
    let zlab_last: Array1<f64> = geometry.zlab_last();

    let acc1 = &mut Accelerator::new();
    let acc2 = &mut Accelerator2d::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = geometry.r_of_psi(psi, acc1).unwrap();
    let _: f64 = geometry.psi_of_r(r, acc1).unwrap();
    let _: f64 = geometry.rlab_of_psi(psi, theta, acc2).unwrap();
    let _: f64 = geometry.zlab_of_psi(psi, theta, acc2).unwrap();

    assert!(matches!(
        geometry.r_of_psip(psip, acc1),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.psip_of_r(r, acc1),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.rlab_of_psip(psip, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.zlab_of_psip(psip, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.jacobian_of_psi(psi, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.jacobian_of_psip(psip, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
}

#[test]
fn nc_geometry() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ1d = Interpolation1dType::Cubic;
    let typ2d = Interpolation2dType::Bicubic;
    let builder = NcGeometryBuilder::new(&path, typ1d, typ2d);
    let geometry = dbg!(builder.build().unwrap());

    assert_eq!(geometry.psi_state(), FluxCoordinateState::Good);
    assert_eq!(geometry.psip_state(), FluxCoordinateState::Good);

    let equilibrium_type: ObjectType = geometry.object_type();
    let netcdf_version: semver::Version = geometry.netcdf_version();
    let path: PathBuf = geometry.path();
    let interp1d_type: Interpolation1dType = geometry.interp1d_type();
    let interp2d_type: Interpolation2dType = geometry.interp2d_type();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let zaxis: f64 = geometry.zaxis();
    let rgeo: f64 = geometry.rgeo();
    let rlast: f64 = geometry.rlast();
    let shape: (usize, usize) = geometry.shape();
    let psi_state: FluxCoordinateState = geometry.psi_state();
    let psip_state: FluxCoordinateState = geometry.psip_state();
    let psi_last: f64 = geometry.psi_last().unwrap();
    let psip_last: f64 = geometry.psip_last().unwrap();
    let psi_array: Array1<f64> = geometry.psi_array().unwrap();
    let psip_array: Array1<f64> = geometry.psip_array().unwrap();
    let theta_array: Array1<f64> = geometry.theta_array();
    let r_array: Array1<f64> = geometry.r_array();
    let rlab_array: Array2<f64> = geometry.rlab_array();
    let zlab_array: Array2<f64> = geometry.rlab_array();
    let jacobian_array: Array2<f64> = geometry.rlab_array();
    let rlab_last: Array1<f64> = geometry.rlab_last();
    let zlab_last: Array1<f64> = geometry.zlab_last();

    let acc1 = &mut Accelerator::new();
    let acc2 = &mut Accelerator2d::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = geometry.psip_of_psi(psi, acc1).unwrap();
    let _: f64 = geometry.psi_of_psip(psip, acc1).unwrap();
    let _: f64 = geometry.r_of_psi(psi, acc1).unwrap();
    let _: f64 = geometry.r_of_psip(psip, acc1).unwrap();
    let _: f64 = geometry.psi_of_r(r, acc1).unwrap();
    let _: f64 = geometry.psip_of_r(r, acc1).unwrap();
    let _: f64 = geometry.rlab_of_psi(psi, theta, acc2).unwrap();
    let _: f64 = geometry.rlab_of_psip(psip, theta, acc2).unwrap();
    let _: f64 = geometry.zlab_of_psi(psi, theta, acc2).unwrap();
    let _: f64 = geometry.zlab_of_psip(psip, theta, acc2).unwrap();
    let _: f64 = geometry.jacobian_of_psi(psi, theta, acc2).unwrap();
    let _: f64 = geometry.jacobian_of_psip(psip, theta, acc2).unwrap();
}
