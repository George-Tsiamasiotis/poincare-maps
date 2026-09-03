//! Test particle orbit classification.
//!
//! These tests are accompanied by Python scripts to visualize the `(E, Pζ)` space and the orbits.

use std::f64::consts::PI;

use dexter_machine::*;
use dexter_simulate::*;

const MU: f64 = 6e-5;

fn create_objects() -> (ParabolicQfactor, LarCurrent, LarBfield) {
    let lcfs = LastClosedFluxSurface::Toroidal(0.03);
    let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    (qfactor, current, bfield)
}

/// CuPassing-Confined inside the left wall parabola.
#[test]
#[rustfmt::skip]
fn orbit_alpha() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.001);
    let pzeta0 = - 0.8 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Alpha);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CuPassing-Lost inside the right wall parabola, outside the left wall parabola and left from the
/// `Pζ/ψp_last = -1` line.
#[test]
#[rustfmt::skip]
fn orbit_beta() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.02);
    let pzeta0 = - 1.5 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Beta);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingLost);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// Trapped-Lost inside the right wall parabola.
#[test]
#[rustfmt::skip]
fn orbit_gamma() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.01);
    let pzeta0 = - 0.8 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Gamma);
    assert_eq!(particle.orbit_type(), OrbitType::TrappedLost);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// CoPassing-Lost inside the right wall parabola, outside the left wall parabola, outside of the
/// axis parabola and outside the trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_delta() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.018);
    let pzeta0 = - 0.6 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, PI, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Delta);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingLost);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// CuPassing-Confined inside the right wall parabola, outside the left wall parabola, outside of the
/// axis parabola and outside the trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_epsilon() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.003);
    let pzeta0 = - 0.6 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, PI, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Epsilon);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CoPassing-Lost inside the right wall parabola and inside axis parabola.
#[test]
#[rustfmt::skip]
fn orbit_zeta() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.025);
    let pzeta0 = - 0.4 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, PI, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Zeta);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingLost);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// CoPassingConfined inside the magnetic axis parabola.
#[test]
#[rustfmt::skip]
fn orbit_eta() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.015);
    let pzeta0 = - 0.1 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Eta);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingConfined);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Potato
#[test]
#[rustfmt::skip]
fn orbit_theta() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.0045);
    let pzeta0 = - 0.0448 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Theta);
    assert_eq!(particle.orbit_type(), OrbitType::Potato);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Trapped-Confined outside of all parabolas.
#[test]
#[rustfmt::skip]
fn orbit_iota() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.016);
    let pzeta0 = - 0.6 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Iota);
    assert_eq!(particle.orbit_type(), OrbitType::TrappedConfined);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CoPassingConfined outside all parabolas and above trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_kappa() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.025);
    let pzeta0 = -0.36 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 0.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Kappa);
    assert_eq!(particle.orbit_type(), OrbitType::CoPassingConfined);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CuPassingConfined outside all parabolas and above trapped-passing boundary.
#[test]
#[rustfmt::skip]
fn orbit_lambda() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.001);
    let pzeta0 = -0.36 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 0.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Lambda);
    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Stagnated
#[test]
#[rustfmt::skip]
fn orbit_mu() {
    let (qfactor, current, bfield )= create_objects();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let psi0 = InitialFlux::Toroidal(0.0014);
    let pzeta0 = - 0.0 * machine.qfactor().psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, MU);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(machine);

    assert_eq!(particle.energy_pzeta_position(), EnergyPzetaPosition::Mu);
    assert_eq!(particle.orbit_type(), OrbitType::Stagnated);

    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}
