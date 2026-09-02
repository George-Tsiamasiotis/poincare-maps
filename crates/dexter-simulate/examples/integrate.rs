//! `Particle::integrate` routine.

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use dexter_equilibrium::{Interpolation1dType::Akima, Interpolation2dType::Bicubic};
use dexter_simulate::*;
use std::path::Path;

fn main() {
    analytical_equilibrium_integration();
    numerical_equilibrium_integration();
}

fn analytical_equilibrium_integration() {
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
    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 7e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    particle.integrate(&equilibrium, (0.0, 1e5), &SolverParams::default());
    particle.print_caches();
    dbg!(&particle);
}

fn numerical_equilibrium_integration() {
    // Equilibrium setup
    let path = Path::new("crates/dexter-simulate").join(TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(vec![
            Box::new(
                NcFluteModeBuilder::new(&path, Akima, 2, 1)
                    .with_phase_method(PhaseMethod::Interpolation)
                    .build()
                    .unwrap(),
            ),
            Box::new(
                NcFluteModeBuilder::new(&path, Akima, 3, 2)
                    .with_phase_method(PhaseMethod::Interpolation)
                    .build()
                    .unwrap(),
            ),
        ]),
    };

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Poloidal(0.2), 0.0, 0.0, 1e-4, 7e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    particle.integrate(&&equilibrium, (0.0, 1e5), &SolverParams::default());
    particle.print_caches();
    dbg!(&particle);
}
