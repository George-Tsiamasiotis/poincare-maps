//! Test Perturbation's functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;

#[test]
#[rustfmt::skip]
fn empty_perturbation() {
    let p = Perturbation::zero();

    let mut caches: DynModeCaches = p.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    assert_eq!(p.p_of_psi(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.p_of_psip(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_dpsi(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_dpsip(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dtheta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dtheta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dzeta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dzeta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dt(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dt(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
}

#[test]
#[rustfmt::skip]
fn cos_toroidal_lcfs_perturbation() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let per = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 1, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    ]);
    assert_eq!(per.count(), 3);

    let mut c: DynModeCaches = per.generate_caches();


    let (p, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);

    let _: f64 = per.p_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_dpsi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psi_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psi_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psi_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(per.p_of_psip(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_dpsip(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psip_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psip_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psip_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero
}

#[test]
#[rustfmt::skip]
fn cos_poloidal_lcfs_perturbation() {
    let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    let per = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 1, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    ]);
    assert_eq!(per.count(), 3);

    let mut c: DynModeCaches = per.generate_caches();


    let (p, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);

    let _: f64 = per.p_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_dpsip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psip_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psip_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psip_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(per.p_of_psi(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_dpsi(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psi_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psi_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psi_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero
}

#[test]
#[rustfmt::skip]
fn nc_perturbation() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Steffen;
    let m1 = NcFluteModeBuilder::new(&path, typ, 2, 1)
        .build()
        .unwrap();
    let m2 = NcFluteModeBuilder::new(&path, typ, 2, 2)
        .build()
        .unwrap();
    let m3 = NcFluteModeBuilder::new(&path, typ, 3, 2)
        .build()
        .unwrap();

    let per = dbg!(Perturbation::new(vec![Box::new(m1), Box::new(m2), Box::new(m3)]));
    assert_eq!(per.count(), 3);

    let mut caches: DynModeCaches = per.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    let _: f64 = per.p_of_psi(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.p_of_psip(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_dpsi(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_dpsip(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psi_dtheta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psip_dtheta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psi_dzeta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psip_dzeta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psi_dt(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psip_dt(psi, theta, zeta, t, &mut caches).unwrap();
}

#[test]
#[rustfmt::skip]
fn mixed_perturbation() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let typ = Interpolation1dType::Cubic;
    let nc_mode = NcFluteModeBuilder::new(&path, typ, 2, 1)
        .build()
        .unwrap();
    let mode = FluteMode::new(1e-3, lcfs, 1, 3, 0.0);

    let per = dbg!(Perturbation::new(vec![Box::new(nc_mode), Box::new(mode)]));
    assert_eq!(per.count(), 2);

    let mut caches: DynModeCaches = per.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    let _: f64 = per.p_of_psi(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_dpsi(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psi_dtheta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psi_dzeta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = per.dp_of_psi_dt(psi, theta, zeta, t, &mut caches).unwrap();
}
