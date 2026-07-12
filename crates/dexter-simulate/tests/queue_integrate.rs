//! Test `Queue::integrate` routine.

#![allow(non_snake_case)]

use dexter_equilibrium::*;
use dexter_simulate::*;
use ndarray::Array1;

#[test]
fn queue_integrate_parQ_larC_larB_cosP() -> Result<(), SimulationError> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(ParabolicQfactor::new(1.1, 1.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(&[
            Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
        ]),
    };

    // Initial Conditions setup
    let particle_count = 10;
    let psis = equilibrium.qfactor.psi_last() * Array1::linspace(0.0, 0.5, particle_count);
    let psis = toroidal_fluxes(&psis.to_vec());
    let initial_conditions = QueueInitialConditions::boozer(
        &vec![0.0; particle_count],
        &psis.to_vec(),
        &vec![0.0; particle_count],
        &vec![0.0; particle_count],
        &vec![1e-6; particle_count],
        &vec![2e-6; particle_count],
    )?;

    let mut queue = Queue::new(&initial_conditions);

    assert_eq!(queue.particle_count(), particle_count);
    assert_eq!(queue.routine(), Routine::None);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Initialized)
            .count(),
        particle_count
    );

    queue.integrate(&equilibrium, (0.0, 1e4), &SolverParams::default());

    // All but the first at `ψ=0` should be integrated.
    assert!(queue[0].integration_status() == IntegrationStatus::OutOfBoundsInitialization);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Integrated)
            .count(),
        particle_count - 1
    );

    println!("{queue:#?}");
    Ok(())
}
