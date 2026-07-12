//! Test particle orbit classification.
//!
//! These tests are accompanied by Python scripts to visualize the `(E, Pζ)` space and the orbits.

use std::f64::consts::PI;

use dexter_equilibrium::*;
use dexter_simulate::*;

const MU: f64 = 6e-5;

fn create_equilibrium() -> Equilibrium {
    let lcfs = LastClosedFluxSurface::Toroidal(0.03);
    Equilibrium {
        geometry: None,
        qfactor: Box::new(ParabolicQfactor::new(1.1, 3.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::zero(),
    }
}

/// CuPassing-Confined inside the left wall parabola.
#[test]
#[rustfmt::skip]
fn orbit_alpha() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.001);
    let pzeta0 = - 0.8 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Alpha);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CuPassing-Lost inside the right wall parabola, outside the left wall parabola and left from the
/// `Pζ/ψp_last = -1` line.
#[test]
#[rustfmt::skip]
fn orbit_beta() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.02);
    let pzeta0 = - 1.5 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Beta);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingLost);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// Trapped-Lost inside the right wall parabola.
#[test]
#[rustfmt::skip]
fn orbit_gamma() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.01);
    let pzeta0 = - 0.8 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Gamma);
    assert_eq!(particle.orbit_type(), OrbitType::TrappedLost);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// CoPassing-Lost inside the right wall parabola, outside the left wall parabola, outside of the
/// axis parabola and outside the trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_delta() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.018);
    let pzeta0 = - 0.6 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, PI, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Delta);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingLost);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// CuPassing-Confined inside the right wall parabola, outside the left wall parabola, outside of the
/// axis parabola and outside the trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_epsilon() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.003);
    let pzeta0 = - 0.6 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, PI, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Epsilon);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CoPassing-Lost inside the right wall parabola and inside axis parabola.
#[test]
#[rustfmt::skip]
fn orbit_zeta() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.025);
    let pzeta0 = - 0.4 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, PI, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Zeta);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingLost);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// CoPassingConfined inside the magnetic axis parabola.
#[test]
#[rustfmt::skip]
fn orbit_eta() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.015);
    let pzeta0 = - 0.1 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Eta);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingConfined);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Potato
#[test]
#[rustfmt::skip]
fn orbit_theta() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.0045);
    let pzeta0 = - 0.0448 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Theta);
    assert_eq!(particle.orbit_type(), OrbitType::Potato);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Trapped-Confined outside of all parabolas.
#[test]
#[rustfmt::skip]
fn orbit_iota() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.016);
    let pzeta0 = - 0.6 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Iota);
    assert_eq!(particle.orbit_type(), OrbitType::TrappedConfined);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CoPassingConfined outside all parabolas and above trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_kappa() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.025);
    let pzeta0 = -0.36 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 0.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Kappa);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingConfined);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CuPassingConfined outside all parabolas and above trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_lambda() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.001);
    let pzeta0 = -0.36 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 0.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Lambda);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Stagnated
#[test]
#[rustfmt::skip]
fn orbit_mu() {
    let equilibrium = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.0014);
    let pzeta0 = - 0.0 * equilibrium.qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&equilibrium);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Mu);
    assert_eq!(particle.orbit_type(), OrbitType::Stagnated);

    particle.close(&equilibrium, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}
