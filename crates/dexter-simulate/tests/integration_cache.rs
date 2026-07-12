//! Tests that the caches work exactly as expected.

use std::path::PathBuf;

use dexter_equilibrium::{extract::TOROIDAL_TEST_NETCDF_PATH, *};
use dexter_simulate::*;

#[test]
fn integration_cache_analytical_eq_flute_mode() {
    use InitialFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(ParabolicQfactor::new(1.1, 3.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(&[
            Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 2, 1, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 2, 2, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 2, 3, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 2, 4, 0.0)),
        ]),
    };

    let initial = InitialConditions::boozer(0.0, Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);

    particle.integrate(&equilibrium, (0.0, 1e5), &SolverParams::default());
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Integrated
    ));
    let steps = particle.steps_taken();
    assert!(steps > 1000);
    let stats = particle.cache_stats().clone();

    // Per mode:
    //     6 evaluations in rkf45
    //          `h`, `dh_dflux` `dh_dtheta`, `dh_dzeta` -> 1 miss, 3 hits
    //
    // Also add 1 miss and 3 hits to account for the state created on setup.

    assert_eq!(stats.mode_cache_misses, 7 * (steps * 6 * 1 + 1));
    assert_eq!(stats.mode_cache_hits, 7 * (steps * 6 * 3 + 3));
    //
    // Each flute mode results in 2 hits and 1 miss for each evaluation
    assert_eq!(stats.mode_cache_hits, 3 * stats.mode_cache_misses);
}

#[test]
fn integration_cache_nc_eq_nc_flute_mode() {
    use InitialFlux::*;
    let path = PathBuf::from(TOROIDAL_TEST_NETCDF_PATH);
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

    let initial = InitialConditions::boozer(0.0, Toroidal(0.2), 0.0, 0.0, 1e-8, 1e-10);
    let mut particle = Particle::new(&initial);

    particle.integrate(
        &equilibrium,
        (0.0, 1e7),
        &SolverParams {
            method: SteppingMethod::FixedStep(2000.0),
            ..Default::default()
        },
    );
    particle.print_cache_stats();
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Integrated
    ));
    let steps = particle.steps_taken();
    assert!(steps > 1000);
    let stats = particle.cache_stats().clone();

    // `ψ` Accelerator
    // We chose a very low energetic particle, so it will follow the same flux surface, so
    // after the first GCState evaluation we have:
    //  2 modes
    //     6 evaluations in rkf45
    //          `other_flux`, `q`, `g`, `i`, `dg`, `di`, `b`, `db_dflux`, `db_dtheta`-> 0 misses 9 hits
    //
    // Also add 1 miss and 8+9 hits (2nd mode is a hit) to account for the state created on setup.
    assert_eq!(stats.psi_acc.misses(), 1);
    assert_eq!(stats.psi_acc.hits(), (8 + 9) + 2 * steps * 6 * 9);

    assert_eq!(stats.psip_acc.hits(), 0);
    assert_eq!(stats.psip_acc.misses(), 0);

    // `θ` Accelerator
    //      2 modes
    //          6 evaluations in rkf45
    //              `b`, `db_dflux`, `db_dtheta`
    //
    // We cannot now how many hits and misses, however their sum must be equal to the total
    // amount of evaluations
    //
    // Also add 6 evaluations to account for the state created on setup.
    //
    assert_eq!(
        stats.theta_acc.hits() + stats.theta_acc.misses(),
        6 + 2 * steps * 6 * 3
    );

    // Per mode:
    //     6 evaluations in rkf45
    //          `h`, `dh_dflux`, `dh_dtheta`, `dh_dzeta` -> 1 miss, 3 hits
    //
    // Also add 1 miss and 3 hits to account for the state created on setup.
    assert_eq!(stats.mode_cache_misses, 2 * (steps * 6 * 1 + 1));
    assert_eq!(stats.mode_cache_hits, 2 * (steps * 6 * 3 + 3));
    assert_eq!(stats.mode_cache_hits, 3 * stats.mode_cache_misses);
}
