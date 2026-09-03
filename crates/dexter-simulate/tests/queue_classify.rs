//! Test `Queue::integrate` routine.

#![allow(non_snake_case)]
#![allow(unused_variables)]

use dexter_machine::*;
use dexter_simulate::*;
use ndarray::Array1;

#[test]
fn queue_classify_parQ_larC_larB_cosP() -> Result<(), SimulationError> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = ParabolicQfactor::new(1.1, 1.9, lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    // Initial Conditions setup
    let particle_count = 10;
    let psis = machine.qfactor().psi_last() * Array1::linspace(0.001, 0.5, particle_count);
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

    queue.classify(machine);

    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Initialized)
            .count(),
        particle_count
    );

    println!("{queue:#?}");

    Ok(())
}
