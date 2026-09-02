//! `Queue::integrate` routine.

use dexter_equilibrium::*;
use dexter_simulate::*;
use ndarray::Array1;

fn main() -> Result<(), SimulationError> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(ParabolicQfactor::new(1.1, 3.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(vec![
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]),
    };

    // Initial Conditions setup
    let particle_count = 200;
    let psis = equilibrium.qfactor.psi_last() * Array1::linspace(0.1, 0.9, particle_count);
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
    queue.integrate(&equilibrium, (0.0, 1e6), &SolverParams::default());
    println!("{queue:#?}");
    Ok(())
}
