//! Test Modes functionality.

#![allow(unused_variables)]

use dexter_machine::*;
use ndarray::Array1;
use std::f64::consts::PI;
use std::path::PathBuf;

#[test]
#[rustfmt::skip]
fn flute_mode_toroidal_lcfs() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let mode = dbg!(FluteMode::new(1e-3, lcfs, 3, 2, PI));
    assert_eq!(mode.psi_state(), FluxCoordinateState::Good);
    assert_eq!(mode.psip_state(), FluxCoordinateState::Bad);

    assert_eq!(mode.machine_type(), MachineType::Analytical);
    assert_eq!(mode.epsilon(), 1e-3);
    assert_eq!(mode.m(), 3);
    assert_eq!(mode.n(), 2);
    assert_eq!(mode.phase(), PI);

    let p = 0.01;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = mode.generate_cache();

    let _: f64 = mode.ampl_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.phase_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.m_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_dpsi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psi_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psi_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psi_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(mode.ampl_of_psip(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.phase_of_psip(p, theta, zeta, t, &mut c).is_ok()); // returns a constant
    assert!(mode.m_of_psip(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_dpsip(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_of_psip_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_of_psip_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_of_psip_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero

    assert_eq!(c.misses(), 1);
    assert_eq!(c.hits(), 4);
}

#[test]
#[rustfmt::skip]
fn flute_mode_poloidal_lcfs() {
    let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    let mode = dbg!(FluteMode::new(1e-3, lcfs, 3, 2, PI));
    assert_eq!(mode.psi_state(), FluxCoordinateState::Bad);
    assert_eq!(mode.psip_state(), FluxCoordinateState::Good);

    assert_eq!(mode.machine_type(), MachineType::Analytical);
    assert_eq!(mode.epsilon(), 1e-3);
    assert_eq!(mode.m(), 3);
    assert_eq!(mode.n(), 2);
    assert_eq!(mode.phase(), PI);

    let p = 0.01;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = mode.generate_cache();

    let _: f64 = mode.ampl_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.phase_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.m_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_dpsip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psip_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psip_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psip_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(mode.ampl_of_psi(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.phase_of_psi(p, theta, zeta, t, &mut c).is_ok()); // returns a constant
    assert!(mode.m_of_psi(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_dpsi(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_of_psi_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_of_psi_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(mode.dm_of_psi_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero

    assert_eq!(c.misses(), 1);
    assert_eq!(c.hits(), 4);
}

#[test]
fn nc_flute_mode() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Steffen;
    let (m, n) = (3, 2);
    use PhaseMethod::Interpolation;
    let builder = NcFluteModeBuilder::new(&path, typ, m, n).with_phase_method(Interpolation);
    let mode = dbg!(builder.build().unwrap());
    assert_eq!(mode.psi_state(), FluxCoordinateState::Good);
    assert_eq!(mode.psip_state(), FluxCoordinateState::Good);

    assert_eq!(mode.machine_type(), MachineType::Numerical);
    assert_eq!(mode.m(), 3);
    assert_eq!(mode.n(), 2);
    assert!(matches!(mode.phase_method(), Interpolation));

    let netcdf_version: semver::Version = mode.netcdf_version();
    let path: PathBuf = mode.path();
    let interp_type: Interpolation1dType = mode.interp_type();
    let psi_state: FluxCoordinateState = mode.psi_state();
    let psip_state: FluxCoordinateState = mode.psip_state();
    let psi_last: f64 = mode.psi_last().unwrap();
    let psip_last: f64 = mode.psip_last().unwrap();
    let psi_array: Array1<f64> = mode.psi_array().unwrap();
    let psip_array: Array1<f64> = mode.psip_array().unwrap();
    let alpha_array: Array1<f64> = mode.alpha_array();
    let phase_array: Array1<f64> = mode.phase_array();

    let phase_average: Option<f64> = mode.phase_average();

    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = mode.generate_cache();

    let _: f64 = mode.ampl_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.ampl_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.phase_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.phase_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.m_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.m_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psi_dtheta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode
        .dm_of_psip_dtheta(psip, theta, zeta, t, &mut c)
        .unwrap();
    let _: f64 = mode.dm_of_psi_dzeta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.dm_of_psip_dzeta(psip, theta, zeta, t, &mut c).unwrap();
    assert_eq!(mode.dm_of_psi_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(
        mode.dm_of_psip_dt(psi, theta, zeta, t, &mut c).unwrap(),
        0.0
    );
}
