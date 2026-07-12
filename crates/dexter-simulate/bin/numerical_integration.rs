//! Integration of a particle for a long time, useful for profiling.

#![expect(clippy::unwrap_used, reason = "not important")]

use dexter_equilibrium::extract::TOROIDAL_TEST_NETCDF_PATH;
use dexter_equilibrium::{
    Equilibrium, NcBfieldBuilder, NcCurrentBuilder, NcFluteModeBuilder, NcQfactorBuilder,
    Perturbation, PhaseMethod,
};
use dexter_simulate::{InitialConditions, InitialFlux, IntegrationStatus, Particle, SolverParams};
use std::path::Path;

fn main() {
    // Equilibrium setup
    let path = Path::new("crates/dexter-simulate").join(TOROIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
    let perturbation = Perturbation::new(&[
        Box::new(
            NcFluteModeBuilder::new(&path, "steffen", 2, 1)
                .with_phase_method(PhaseMethod::Interpolation)
                .build()
                .unwrap(),
        ),
        Box::new(
            NcFluteModeBuilder::new(&path, "steffen", 3, 2)
                .with_phase_method(PhaseMethod::Interpolation)
                .build()
                .unwrap(),
        ),
    ]);

    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(qfactor),
        bfield: Box::new(bfield),
        current: Box::new(current),
        perturbation,
    };

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 1.0, 0.0, 1e-4, 1e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    let teval = (0.0, 1e10);
    particle.integrate(&equilibrium, teval, &SolverParams::default());
    dbg!(&particle);
    particle.print_cache_stats();
    assert!(
        matches!(
            particle.integration_status(),
            IntegrationStatus::TimedOut(..)
        ),
        "particle is supposed to time out"
    );
}
