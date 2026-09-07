//! `θ-ψ` orbit closing and qkinetic calculation checks.

#![allow(non_snake_case)]

mod common;

use approx::assert_relative_eq as ar;
use dexter_machine::*;
use dexter_simulate::*;
use std::f64::consts::TAU;

/// Simplest case
/// Use a very low `ρ` and `μ=0` to track a magnetic field line, which should have
/// `qkinetic = qmagnetic = 1`
///
/// Note that this orbit has essentially Δψ=0. The high order hamiltonian terms are essentially
/// zero, so the `first` and `final` values should be almost equal, up to the solver's tolerance.
///
/// Results calculated manually with the `integrate()` routine.
#[test]
#[rustfmt::skip]
fn field_line_single_period_uniQ() {
    let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.5));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.3), 1.0, 0.0, 1e-12, 0.0);
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &SolverParams::default());
    assert!(matches!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1)));

    let expected_dt = 17085113170546.295;
    let expected_omega_theta = TAU / expected_dt;
    let expected_omega_zeta = TAU / expected_dt;

    let eps=1e-10;
    ar!(particle.t_array().last().copied().unwrap(), expected_dt, epsilon=eps);
    ar!(particle.psi_array().first().copied().unwrap(), particle.psi_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.psip_array().first().copied().unwrap(), particle.psip_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.theta_array().first().copied().unwrap(), particle.theta_array().last().copied().unwrap() - TAU, epsilon=eps);
    ar!(particle.zeta_array().first().copied().unwrap(), particle.zeta_array().last().copied().unwrap() - TAU, epsilon=eps);
    ar!(particle.rho_array().first().copied().unwrap(), particle.rho_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.pzeta_array().first().copied().unwrap(), particle.pzeta_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.ptheta_array().first().copied().unwrap(), particle.ptheta_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.frequencies().omega_theta.unwrap(), expected_omega_theta, epsilon=eps);
    ar!(particle.frequencies().omega_zeta.unwrap(), expected_omega_zeta, epsilon=eps);
    ar!(particle.frequencies().qkinetic.unwrap(), 1.0, epsilon=eps);
}

#[test]
#[rustfmt::skip]
fn trapped_particle_single_period_uniQ() {
    let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.02), 1.0, 0.0, 1e-6, 1e-6);
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &SolverParams::default());
    particle.classify(machine);
    assert!(matches!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1)));
    assert_eq!(particle.orbit_type(), OrbitType::TrappedConfined);

    let expected_dt = 17708.921170693902;
    let expected_dzeta = 0.07681124462891854; // from manual integration
    let expected_omega_theta = TAU / expected_dt;
    let expected_omega_zeta = expected_dzeta / expected_dt;

    let eps = 1e-6;
    ar!(particle.t_array().last().copied().unwrap(), expected_dt, epsilon=eps);
    ar!(particle.psi_array().first().copied().unwrap(), particle.psi_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.theta_array().first().copied().unwrap(), particle.theta_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.ptheta_array().first().copied().unwrap(), particle.ptheta_array().last().copied().unwrap(), epsilon=eps);
    ar!(particle.frequencies().omega_theta.unwrap(), expected_omega_theta, epsilon=eps);
    ar!(particle.frequencies().omega_zeta.unwrap(), expected_omega_zeta, epsilon=eps);
    ar!(particle.frequencies().qkinetic.unwrap(), 0.012224889267733182, epsilon=eps); // derived
}
