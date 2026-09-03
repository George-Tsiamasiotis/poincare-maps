//! Integration of a particle for a long time, useful for profiling.

use dexter_machine::{
    FluteMode, LarBfield, LarCurrent, LastClosedFluxSurface, MachineBuilder, ParabolicQfactor,
    Perturbation,
};
use dexter_simulate::{
    InitialConditions, InitialFlux, IntegrationStatus, Particle, SolverParams, SteppingMethod,
};

fn main() {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    let teval = (0.0, 10.0);
    let solver_params = SolverParams {
        method: SteppingMethod::EnergyAdaptiveStep,
        energy_rel_tol: 1e-17,
        energy_abs_tol: 1e-18,
        max_steps: 10_000_000,
        ..Default::default()
    };
    particle.integrate(machine, teval, &solver_params);
    dbg!(&particle);
    assert!(
        matches!(
            particle.integration_status(),
            IntegrationStatus::TimedOut(..)
        ),
        "particle is supposed to time out"
    );
}
