//! Different integration `SteppingMethods` testing.

mod common;

use crate::common::check_integrated_particle_arrays;
use dexter_machine::*;
use dexter_simulate::*;
use ndarray::Axis;

#[test]
#[rustfmt::skip]
fn different_stepping_methods() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.5));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).with_perturbation(&perturbation).build();

    let energy_params =  SolverParams{
            method: SteppingMethod::EnergyAdaptiveStep,
            energy_rel_tol: 1e-13,
            energy_abs_tol: 1e-15,
            ..Default::default()
        };
    let error_params =  SolverParams{
            method: SteppingMethod::ErrorAdaptiveStep,
            error_rel_tol: 1e-13,
            energy_abs_tol: 1e-15,
            ..Default::default()
        };

    let fixed_params =  SolverParams{
            method: SteppingMethod::FixedStep(10.0),
            ..Default::default()
        };

    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.3), 0.0, 0.0, 1e-4, 1e-6);

    let mut energy_particle = Particle::new(&initial);
    let mut error_particle = Particle::new(&initial);
    let mut fixed_particle = Particle::new(&initial);

    let teval = (0.0, 1e5);
    energy_particle.integrate(machine, teval, &energy_params);
    error_particle.integrate(machine, teval, &error_params);
    fixed_particle.integrate(machine, teval, &fixed_params);

    dbg!(&energy_particle);
    dbg!(&error_particle);
    dbg!(&fixed_particle);

    assert!(matches!(energy_particle.integration_status(), IntegrationStatus::Integrated));
    assert!(matches!(error_particle.integration_status(), IntegrationStatus::Integrated));
    assert!(matches!(fixed_particle.integration_status(), IntegrationStatus::Integrated));

    assert!(fixed_particle.t_array().diff(1, Axis(0)).iter().all(|d| *d == 10.0));

    check_integrated_particle_arrays(&energy_particle);
    check_integrated_particle_arrays(&error_particle);
    check_integrated_particle_arrays(&fixed_particle);
}
