//! `Particle::intersect` routine.

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use dexter_simulate::*;
use std::path::Path;

fn main() {
    analytical_equilibrium_intersect();
    numerical_equilibrium_intersect();
}

fn analytical_equilibrium_intersect() {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(ParabolicQfactor::new(1.1, 3.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(&[
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]),
    };

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 0.0, 100);
    let mut particle = Particle::new(&initial);

    // Calculate intersections
    particle.intersect(&equilibrium, &intersect_params, &SolverParams::default());
    dbg!(&particle);
}

fn numerical_equilibrium_intersect() {
    // Equilibrium setup
    let path = Path::new("crates/dexter-simulate").join(TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, "steffen").build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, "steffen").build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, "bicubic").build().unwrap()),
        perturbation: Perturbation::new(&[
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
        ]),
    };

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-6, 0.0);
    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 0.0, 100);
    let mut particle = Particle::new(&initial);

    // Calculate intersections
    particle.intersect(&equilibrium, &intersect_params, &SolverParams::default());
    dbg!(&particle);
}
