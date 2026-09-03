//! Comparison of the different integration `SteppingMethods`.

use dexter_machine::*;
use dexter_simulate::*;

fn main() {
    use InitialFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = UnityQfactor::new(lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-2, lcfs, 1, 1, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let energy_params = SolverParams {
        method: SteppingMethod::EnergyAdaptiveStep,
        energy_rel_tol: 1e-14,
        energy_abs_tol: 1e-15,
        ..Default::default()
    };
    let error_params = SolverParams {
        method: SteppingMethod::ErrorAdaptiveStep,
        error_rel_tol: 1e-22,
        error_abs_tol: 1e-23,
        ..Default::default()
    };

    let fixed_params = SolverParams {
        method: SteppingMethod::FixedStep(4.0),
        ..Default::default()
    };

    let initial = InitialConditions::boozer(0.0, Toroidal(0.3), 0.0, 0.0, 1e-4, 1e-6);

    let mut energy_particle = Particle::new(&initial);
    let mut error_particle = Particle::new(&initial);
    let mut fixed_particle = Particle::new(&initial);

    let teval = (0.0, 1e5);
    energy_particle.integrate(machine, teval, &energy_params);
    error_particle.integrate(machine, teval, &error_params);
    fixed_particle.integrate(machine, teval, &fixed_params);

    println!("Energy adaptive step:");
    print_results(&energy_particle);
    println!("Error adaptive step:");
    print_results(&error_particle);
    println!("Fixed step:");
    print_results(&fixed_particle);
}

fn print_results(particle: &Particle) {
    println!("\tSteps taken: {}", particle.steps_taken());
    println!("\tEnergy variance: {:?}", particle.energy_var().unwrap());
    println!("\tDuration: {:?}", particle.duration());
}
