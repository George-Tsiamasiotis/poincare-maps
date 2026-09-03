//! Compare `Queue::classify` and `Queue::classify_common_mu` routines

use std::time::Instant;

use dexter_machine::*;
use dexter_simulate::*;
use ndarray::Array1;

fn main() -> Result<(), SimulationError> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    // Initial Conditions setup
    let particle_count = 1_000_000;
    let pzetas = machine.qfactor().psip_last() * Array1::linspace(-1.4, 0.2, particle_count);
    let psis = machine.qfactor().psi_last() * Array1::linspace(0.001, 0.5, particle_count);
    let psis = toroidal_fluxes(&psis.to_vec());
    let initial_conditions = QueueInitialConditions::mixed(
        Array1::zeros(particle_count).as_slice().unwrap(),
        &psis.to_vec(),
        Array1::ones(particle_count).as_slice().unwrap(),
        Array1::zeros(particle_count).as_slice().unwrap(),
        pzetas.as_slice().unwrap(),
        Array1::from_elem(particle_count, 7e-6).as_slice().unwrap(),
    )?;

    // ===========================================================================================

    let mut no_common_queue = Queue::new(&initial_conditions);
    let no_common_start = Instant::now();
    no_common_queue.classify(machine);
    let no_common_elapsed = no_common_start.elapsed();

    let mut common_queue = Queue::new(&initial_conditions);
    let common_start = Instant::now();
    common_queue.classify_common_mu(machine);
    let common_elapsed = common_start.elapsed();

    // Sanity check
    for particle_index in 0..no_common_queue.particle_count() {
        assert_eq!(
            no_common_queue[particle_index].orbit_type(),
            common_queue[particle_index].orbit_type(),
            "orbit types must be the same"
        );
    }

    println!("===============================================================");
    println!("Number of particles: {particle_count}");
    println!("`Queue::classify` duration: {no_common_elapsed:?}");
    println!("`Queue::classify_common_mu` duration: {common_elapsed:?}");
    println!("===============================================================");

    Ok(())
}
