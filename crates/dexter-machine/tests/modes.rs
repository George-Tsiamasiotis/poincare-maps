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

    let psi = MagneticFlux::Toroidal(0.01);
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = mode.generate_cache();

    let _: f64 = mode.eval_amplitude(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_phase(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_m(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_flux(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_theta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_zeta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_t(psi, theta, zeta, t, &mut c).unwrap();

    let psip = MagneticFlux::Poloidal(0.01);
    assert!(mode.eval_amplitude(psip, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_phase(psip, theta, zeta, t, &mut c).is_ok()); // returns a constant
    assert!(mode.eval_m(psip, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_flux(psip, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_theta(psip, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_zeta(psip, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_t(psip, theta, zeta, t, &mut c).is_ok()); // returns zero

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

    let psip = MagneticFlux::Poloidal(0.01);
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = mode.generate_cache();

    let _: f64 = mode.eval_amplitude(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_phase(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_m(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_flux(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_theta(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_zeta(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_t(psip, theta, zeta, t, &mut c).unwrap();

    let psi = MagneticFlux::Toroidal(0.01);
    assert!(mode.eval_amplitude(psi, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_phase(psi, theta, zeta, t, &mut c).is_ok()); // returns a constant
    assert!(mode.eval_m(psi, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_flux(psi, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_theta(psi, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_zeta(psi, theta, zeta, t, &mut c).is_err());
    assert!(mode.eval_deriv_t(psi, theta, zeta, t, &mut c).is_ok()); // returns zero

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

    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = mode.generate_cache();

    let _: f64 = mode.eval_amplitude(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_amplitude(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_phase(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_phase(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_m(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_m(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_theta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_theta(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_zeta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = mode.eval_deriv_zeta(psip, theta, zeta, t, &mut c).unwrap();
    assert_eq!(mode.eval_deriv_t(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(mode.eval_deriv_t(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
}
