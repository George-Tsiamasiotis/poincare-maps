//! Guiding Center fixed interval integration.

#![allow(non_snake_case)]

mod common;

use crate::common::check_integrated_particle_arrays;
use approx::{RelativeEq, assert_abs_diff_eq, assert_relative_eq};
use dexter_equilibrium::{
    Interpolation1dType::Akima,
    Interpolation2dType::Bicubic,
    extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH, TOROIDAL_TEST_NETCDF_PATH},
    *,
};
use dexter_simulate::*;
use ndarray::Array1;
use std::{f64::consts::TAU, path::PathBuf};

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_uniQ_larC_larB_cosP() {
    use InitialFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(UnityQfactor::new(lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(&[
            Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
        ]),
    };

    let solver_params = SolverParams::default();

    let initial = InitialConditions::boozer(0.0, Toroidal(0.3), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(&equilibrium, teval, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_ncdQ_ncdC_ncdB_ncdP() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TOROIDAL_TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor : Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current : Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield : Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(&[
            Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, Akima, 3, 2).with_phase_method(Interpolation).build().unwrap()),
        ]),
    };

    let solver_params = SolverParams::default();

    let initial = InitialConditions::boozer(0.0, Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(&equilibrium, teval, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_poloidal_integration_ncdQ_ncdC_ncdB_ncdP() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor : Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current : Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield : Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(&[
            Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, Akima, 3, 2).with_phase_method(Interpolation).build().unwrap()),
        ]),
    };

    let solver_params = SolverParams::default();

    let initial = InitialConditions::boozer(0.0, Poloidal(0.2), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(&equilibrium, teval, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_gcmotion_check_uniQ_larC_larB_cosP() {
    use InitialFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor : Box::new(UnityQfactor::new(lcfs)),
        current : Box::new(LarCurrent::new()),
        bfield : Box::new(LarBfield::new()),
        perturbation: Perturbation::zero(),
    };

    let solver_params = SolverParams::default();

    let psi0 = 0.018365472910927463;
    let initial = InitialConditions::boozer(0.0, Toroidal(psi0), 0.0, 0.0, -0.006634527089072539, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Initialized
    ));

    let teval = (0.0, 3e3);
    particle.integrate(&equilibrium, teval, &solver_params);

    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Integrated
    ));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-18);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    assert_relative_eq!(particle.initial_energy().unwrap(), 1.5189224863170239e-05, epsilon = 1e-20);
    assert_relative_eq!(initial.mu0() * equilibrium.bfield.b_of_psi(psi0, 0.0, &mut Accelerator2d::new()).unwrap(),
        8.083468084746437e-07,
        epsilon = 1e-15
    );
    assert_relative_eq!(*particle.psi_array().last().unwrap(), 0.02059090128879661, epsilon = 1e-5);
    assert_abs_diff_eq!(*particle.theta_array().last().unwrap(), -15.887142479457443, epsilon = 0.1);
    assert_abs_diff_eq!(*particle.zeta_array().last().unwrap(), -15.80106331080741, epsilon = 0.1);
    assert_relative_eq!(*particle.rho_array().last().unwrap(), -0.004409098711254563, epsilon = 1e-5);
    assert_relative_eq!(*particle.ptheta_array().last().unwrap(), 0.020590901288745446, epsilon = 1e-5);
    assert_relative_eq!(particle.pzeta_array().var(1.0), 0.0, epsilon = 1e-10);
    assert_relative_eq!(*particle.pzeta_array().first().unwrap(), -0.025, epsilon = 1e-4);
    assert_relative_eq!(*particle.pzeta_array().last().unwrap(), -0.025, epsilon = 1e-4);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_poloidal_equivalence() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor : Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current : Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield : Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(&[
            Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, Akima, 3, 2).with_phase_method(Interpolation).build().unwrap()),
        ]),
    };

    let solver_params = SolverParams::default();

    let psi0 = 0.2;
    let psip0 = equilibrium.qfactor.psip_of_psi(psi0, &mut Accelerator::new()).unwrap();
    let tor_initial = InitialConditions::boozer(0.0, Toroidal(psi0), 0.0, 0.0, 1e-4, 1e-6);
    let pol_initial = InitialConditions::boozer(0.0, Poloidal(psip0), 0.0, 0.0, 1e-4, 1e-6);

    let mut tor_particle = Particle::new(&tor_initial);
    let mut pol_particle = Particle::new(&pol_initial);
    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Initialized));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 4.14e5);
    tor_particle.integrate(&equilibrium, teval, &solver_params);
    pol_particle.integrate(&equilibrium, teval, &solver_params);
    dbg!(&tor_particle);
    dbg!(&pol_particle);

    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Integrated));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Integrated));

    // Make sure we dont integrate the same particle twice. The energy should have a small
    // arithmetic error due to the ψ-ψp interpolation.
    assert_ne!(tor_particle.initial_energy(), pol_particle.initial_energy());

    assert!(tor_particle.steps_stored() > 10_000);
    assert!(pol_particle.steps_stored() > 10_000);
    assert!(tor_particle.steps_stored() < 10_100);
    assert!(pol_particle.steps_stored() < 10_100);

    check_integrated_particle_arrays(&tor_particle);
    check_integrated_particle_arrays(&pol_particle);

    assert_relative_eq!(tor_particle.initial_energy().unwrap(), pol_particle.initial_energy().unwrap(), epsilon=1e-12);
    assert_relative_eq!(tor_particle.final_energy().unwrap(), pol_particle.final_energy().unwrap(), epsilon=1e-11);

    assert_abs_diff_eq!(
        tor_particle.t_array().last().copied().unwrap(),
        pol_particle.t_array().last().copied().unwrap(),
        epsilon = teval.1 * 1e-3
    );
    assert_abs_diff_eq!(
        tor_particle.psi_array().last().copied().unwrap(),
        pol_particle.psi_array().last().copied().unwrap(),
        epsilon = 1e-3*equilibrium.qfactor.psi_last()
    );
    assert_abs_diff_eq!(
        tor_particle.psip_array().last().copied().unwrap(),
        pol_particle.psip_array().last().copied().unwrap(),
        epsilon = 1e-3*equilibrium.qfactor.psip_last()
    );
    assert_abs_diff_eq!(
        tor_particle.rho_array().last().copied().unwrap(),
        pol_particle.rho_array().last().copied().unwrap(),
        epsilon = 1e-2*tor_initial.rho0().unwrap() // This deviation might be because rho varies quite fast
    );
    assert_abs_diff_eq!(
        tor_particle.theta_array().last().copied().unwrap(),
        pol_particle.theta_array().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_abs_diff_eq!(
        tor_particle.zeta_array().last().copied().unwrap(),
        pol_particle.zeta_array().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_relative_eq!(
        tor_particle.ptheta_array().last().copied().unwrap(),
        pol_particle.ptheta_array().last().copied().unwrap(),
        epsilon = 1e-5
    );
    assert_relative_eq!(
        tor_particle.pzeta_array().last().copied().unwrap(),
        pol_particle.pzeta_array().last().copied().unwrap(),
        epsilon = 1e-5
    );
}

#[test]
#[rustfmt::skip]
fn gc_mixed_boozer_equivalence() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor : Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current : Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield : Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(&[
            Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, Akima, 3, 2).with_phase_method(Interpolation).build().unwrap()),
        ]),
    };

    let solver_params = SolverParams::default();

    let rho0 = 1e-4;
    let psip0 = 0.2;
    let g0 = equilibrium.current.g_of_psip(psip0, &mut Accelerator::new()).unwrap();
    let pzeta0 = rho0 * g0 - psip0;

    let boozer_initial = InitialConditions::boozer(0.0, Poloidal(psip0), 0.0, 0.0, rho0, 1e-6);
    let mixed_initial = InitialConditions::mixed(0.0, Poloidal(psip0), 0.0, 0.0, pzeta0, 1e-6);

    let mut boozer_particle = Particle::new(&boozer_initial);
    let mut mixed_particle = Particle::new(&mixed_initial);
    assert!(matches!(boozer_particle.integration_status(), IntegrationStatus::Initialized));
    assert!(matches!(mixed_particle.integration_status(), IntegrationStatus::PartlyInitialized));

    let teval = (0.0, 2.3e5);
    boozer_particle.integrate(&equilibrium, teval, &solver_params);
    mixed_particle.integrate(&equilibrium, teval, &solver_params);
    dbg!(&boozer_particle);
    dbg!(&mixed_particle);

    assert!(matches!(boozer_particle.integration_status(), IntegrationStatus::Integrated));
    assert!(matches!(mixed_particle.integration_status(), IntegrationStatus::Integrated));

    // Make sure we dont integrate the same particle twice. The energy should have a small
    // arithmetic error due to the ψ-ψp interpolation.
    assert_ne!(boozer_particle.initial_energy(), mixed_particle.initial_energy());

    assert!(boozer_particle.steps_stored() > 10_000);
    assert!(mixed_particle.steps_stored() > 10_000);
    assert!(boozer_particle.steps_stored() < 10_100);
    assert!(mixed_particle.steps_stored() < 10_100);

    check_integrated_particle_arrays(&boozer_particle);
    check_integrated_particle_arrays(&mixed_particle);

    assert_relative_eq!(boozer_particle.initial_energy().unwrap(), mixed_particle.initial_energy().unwrap(), epsilon=1e-12);
    assert_relative_eq!(boozer_particle.final_energy().unwrap(), mixed_particle.final_energy().unwrap(), epsilon=1e-11);

    assert_abs_diff_eq!(
        boozer_particle.t_array().last().copied().unwrap(),
        mixed_particle.t_array().last().copied().unwrap(),
        epsilon = teval.1 * 1e-3
    );
    assert_abs_diff_eq!(
        boozer_particle.psi_array().last().copied().unwrap(),
        mixed_particle.psi_array().last().copied().unwrap(),
        epsilon = 1e-3*equilibrium.qfactor.psi_last()
    );
    assert_abs_diff_eq!(
        boozer_particle.psip_array().last().copied().unwrap(),
        mixed_particle.psip_array().last().copied().unwrap(),
        epsilon = 1e-3*equilibrium.qfactor.psip_last()
    );
    assert_abs_diff_eq!(
        boozer_particle.rho_array().last().copied().unwrap(),
        mixed_particle.rho_array().last().copied().unwrap(),
        epsilon = 1e-2*boozer_initial.rho0().unwrap() // This deviation might be because rho varies quite fast
    );
    assert_abs_diff_eq!(
        boozer_particle.theta_array().last().copied().unwrap(),
        mixed_particle.theta_array().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_abs_diff_eq!(
        boozer_particle.zeta_array().last().copied().unwrap(),
        mixed_particle.zeta_array().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_relative_eq!(
        boozer_particle.ptheta_array().last().copied().unwrap(),
        mixed_particle.ptheta_array().last().copied().unwrap(),
        epsilon = 1e-5
    );
    assert_relative_eq!(
        boozer_particle.pzeta_array().last().copied().unwrap(),
        mixed_particle.pzeta_array().last().copied().unwrap(),
        epsilon = 1e-5
    );
}

#[test]
#[rustfmt::skip]
fn gc_const_pzeta_toroidal_integration_uniQ_larC_larB_cosP() {
    use InitialFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(UnityQfactor::new(lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::zero(),
    };

    let solver_params = SolverParams::default();
    let initial = InitialConditions::mixed(0.0, Toroidal(0.03), 0.0, 0.0, -0.03, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::PartlyInitialized));

    let teval = (0.0, 6e5);
    particle.integrate(&equilibrium, teval, &solver_params);
    let pzeta_array = particle.pzeta_array();
    dbg!(&particle);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);
    assert!(
        pzeta_array.relative_eq(
            &Array1::from_elem(particle.steps_stored(), initial.pzeta0().unwrap()),
            1e-13,
            1e-15
        )
    );
    assert!(pzeta_array.std(1.0) <= 1e-15);

    check_integrated_particle_arrays(&particle);
}
