//! Test `Queue::close` routine.

#![allow(non_snake_case)]
#![allow(unused_variables)]

use dexter_machine::*;
use dexter_simulate::*;
use ndarray::Array1;

#[test]
fn queue_close_parQ_larC_larB_cosP() -> Result<(), SimulationError> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = ParabolicQfactor::new(1.1, 1.9, LastClosedFluxSurface::Toroidal(0.45));
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

    queue.close(machine, 1, &SolverParams::default());

    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::ClosedPeriods(1))
            .count(),
        particle_count
    );

    println!("{queue:#?}");

    let steps_taken = queue.steps_taken_array();
    let steps_stored = queue.steps_stored_array();
    let energy_array = queue.energy_array();
    let qkinetic_array = queue.qkinetic_array();
    let omega_theta_array = queue.omega_theta_array();
    let omega_zeta_array = queue.omega_zeta_array();
    let durations = queue.durations();

    assert!(steps_taken.iter().all(|steps| *steps > 0));
    assert!(steps_stored.iter().all(|steps| *steps == 0)); // close() discards arrays
    assert!(energy_array.iter().all(|energy| energy.is_finite()));
    assert!(qkinetic_array.iter().all(|qkinetic| qkinetic.is_finite()));
    assert!(omega_theta_array.iter().all(|omega| omega.is_finite()));
    assert!(omega_zeta_array.iter().all(|omega| omega.is_finite()));
    assert!(durations.iter().all(|duration| !duration.is_zero()));

    Ok(())
}
