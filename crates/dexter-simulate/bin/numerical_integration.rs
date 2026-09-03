//! Integration of a particle for a long time, useful for profiling.

#![expect(clippy::unwrap_used, reason = "not important")]

use dexter_machine::MachineBuilder;
use dexter_machine::extract::TOROIDAL_TEST_NETCDF_PATH;
use dexter_machine::{
    Interpolation1dType::Steffen, Interpolation2dType::Bicubic, NcBfieldBuilder, NcCurrentBuilder,
    NcFluteModeBuilder, NcQfactorBuilder, Perturbation, PhaseMethod,
};
use dexter_simulate::{InitialConditions, InitialFlux, IntegrationStatus, Particle, SolverParams};
use std::path::Path;

fn main() {
    // Equilibrium setup
    let path = Path::new("crates/dexter-simulate").join(TOROIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, Steffen).build().unwrap();
    let current = NcCurrentBuilder::new(&path, Steffen).build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
    let perturbation = Perturbation::new(vec![
        Box::new(
            NcFluteModeBuilder::new(&path, Steffen, 2, 1)
                .with_phase_method(PhaseMethod::Interpolation)
                .build()
                .unwrap(),
        ),
        Box::new(
            NcFluteModeBuilder::new(&path, Steffen, 3, 2)
                .with_phase_method(PhaseMethod::Interpolation)
                .build()
                .unwrap(),
        ),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 1.0, 0.0, 1e-4, 1e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    let teval = (0.0, 1e10);
    particle.integrate(machine, teval, &SolverParams::default());
    particle.print_caches();
    dbg!(&particle);
    assert!(
        matches!(
            particle.integration_status(),
            IntegrationStatus::TimedOut(..)
        ),
        "particle is supposed to time out"
    );
}
