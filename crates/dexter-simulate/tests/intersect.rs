//! `Particle::intersect` testing.

#![allow(non_snake_case)]

mod common;

use crate::common::check_integrated_particle_arrays;
use approx::*;
use dexter_equilibrium::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
use dexter_equilibrium::*;
use dexter_equilibrium::{Interpolation1dType::Akima, Interpolation2dType::Bicubic};
use dexter_simulate::*;
use std::path::PathBuf;

#[test]
#[rustfmt::skip]
fn gc_toroidal_intersect_uniQ_larC_larB_noP() {
    use InitialFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(UnityQfactor::new(lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::zero(),
    };

    let initial = InitialConditions::boozer(0.0, Toroidal(0.02), 3.14, 0.0, 1e-4, 1e-6);
    let solver_params = SolverParams::default();

    // ConstZeta
    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 3.14, 10);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    particle.intersect(&equilibrium, &intersect_params, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    assert_eq!(particle.steps_stored(), 10);
    assert!(particle.steps_taken() > 10);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    check_integrated_particle_arrays(&particle);

    // ConstTheta
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 3.14, 10);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    particle.intersect(&equilibrium, &intersect_params, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    assert_eq!(particle.steps_stored(), 10);
    assert!(particle.steps_taken() > 10);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_poloidal_intersect_ncdQ_ncdC_ncdB_noP() {
    use InitialFlux::*;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::zero(),
    };

    let initial = InitialConditions::boozer(0.0, Poloidal(0.1), 3.14, 0.0, 1e-4, 1e-6);
    let solver_params = SolverParams::default();

    // ConstZeta
    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 3.14, 3);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    particle.intersect(&equilibrium, &intersect_params, &solver_params);
    dbg!(&particle);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    assert!(particle.steps_stored() == 3);
    assert!(particle.steps_taken() > 3);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    check_integrated_particle_arrays(&particle);

    // ConstTheta
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 3.14, 3);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    particle.intersect(&equilibrium, &intersect_params, &solver_params);
    dbg!(&particle);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    assert!(particle.steps_stored() == 3);
    assert!(particle.steps_taken() > 3);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_poloidal_equivalence_const_theta() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(vec![
            Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, Akima, 3, 2).with_phase_method(Interpolation).build().unwrap()),
        ]),
    };

    let solver_params = SolverParams::default();

    let psi0 = 0.2;
    let psip0 = equilibrium.qfactor.psip_of_psi(psi0, &mut Accelerator::new()).unwrap();
    let tor_initial = InitialConditions::boozer(0.0, Toroidal(psi0), 0.0, 0.0, 1e-4, 0.0);
    let pol_initial = InitialConditions::boozer(0.0, Poloidal(psip0), 0.0, 0.0, 1e-4, 0.0);

    let mut tor_particle = Particle::new(&tor_initial);
    let mut pol_particle = Particle::new(&pol_initial);
    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Initialized));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Initialized));


    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 1.0, 10);
    tor_particle.intersect(&equilibrium, &intersect_params, &solver_params);
    pol_particle.intersect(&equilibrium, &intersect_params, &solver_params);
    dbg!(&tor_particle);
    dbg!(&pol_particle);

    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Intersected));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Intersected));

    assert_eq!(tor_particle.steps_stored(), 10);
    assert_eq!(pol_particle.steps_stored(), 10);

    check_integrated_particle_arrays(&tor_particle);
    check_integrated_particle_arrays(&pol_particle);

    assert_relative_eq!(tor_particle.initial_energy().unwrap(), pol_particle.initial_energy().unwrap(), epsilon=1e-12);
    assert_relative_eq!(tor_particle.final_energy().unwrap(), pol_particle.final_energy().unwrap(), epsilon=1e-11);

    assert!(tor_particle.t_array().relative_eq(&pol_particle.t_array(), 1e-6, 1e-4));
    assert!(tor_particle.psi_array().relative_eq(&pol_particle.psi_array(), 1e-8, 1e-5));
    assert!(tor_particle.psip_array().relative_eq(&pol_particle.psip_array(), 1e-8, 1e-5));
    assert!(tor_particle.theta_array().relative_eq(&pol_particle.theta_array(), 1e-8, 1e-5));
    assert!(tor_particle.zeta_array().relative_eq(&pol_particle.zeta_array(), 1e-8, 1e-5));
    assert!(tor_particle.rho_array().relative_eq(&pol_particle.rho_array(), 1e-8, 1e-5));
    assert!(tor_particle.ptheta_array().relative_eq(&pol_particle.ptheta_array(), 1e-8, 1e-5));
    assert!(tor_particle.pzeta_array().relative_eq(&pol_particle.pzeta_array(), 1e-8, 1e-5));
    assert!(tor_particle.energy_array().relative_eq(&pol_particle.energy_array(), 1e-8, 1e-5));
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_poloidal_equivalence_const_zeta() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let equilibrium = Equilibrium {
        geometry: None,
        qfactor: Box::new(NcQfactorBuilder::new(&path, Akima).build().unwrap()),
        current: Box::new(NcCurrentBuilder::new(&path, Akima).build().unwrap()),
        bfield: Box::new(NcBfieldBuilder::new(&path, Bicubic).build().unwrap()),
        perturbation: Perturbation::new(vec![
            Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
            Box::new(NcFluteModeBuilder::new(&path, Akima, 3, 2).with_phase_method(Interpolation).build().unwrap()),
        ]),
    };

    let solver_params = SolverParams::default();

    let psi0 = 0.2;
    let psip0 = equilibrium.qfactor.psip_of_psi(psi0, &mut Accelerator::new()).unwrap();
    let tor_initial = InitialConditions::boozer(0.0, Toroidal(psi0), 0.0, 0.0, 1e-4, 0.0);
    let pol_initial = InitialConditions::boozer(0.0, Poloidal(psip0), 0.0, 0.0, 1e-4, 0.0);

    let mut tor_particle = Particle::new(&tor_initial);
    let mut pol_particle = Particle::new(&pol_initial);
    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Initialized));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Initialized));


    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 1.0, 10);
    tor_particle.intersect(&equilibrium, &intersect_params, &solver_params);
    pol_particle.intersect(&equilibrium, &intersect_params, &solver_params);
    dbg!(&tor_particle);
    dbg!(&pol_particle);

    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Intersected));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Intersected));

    assert_eq!(tor_particle.steps_stored(), 10);
    assert_eq!(pol_particle.steps_stored(), 10);

    check_integrated_particle_arrays(&tor_particle);
    check_integrated_particle_arrays(&pol_particle);

    assert_relative_eq!(tor_particle.initial_energy().unwrap(), pol_particle.initial_energy().unwrap(), epsilon=1e-12);
    assert_relative_eq!(tor_particle.final_energy().unwrap(), pol_particle.final_energy().unwrap(), epsilon=1e-11);

    assert!(tor_particle.t_array().relative_eq(&pol_particle.t_array(), 1e-6, 1e-4));
    assert!(tor_particle.psi_array().relative_eq(&pol_particle.psi_array(), 1e-8, 1e-5));
    assert!(tor_particle.psip_array().relative_eq(&pol_particle.psip_array(), 1e-8, 1e-5));
    assert!(tor_particle.theta_array().relative_eq(&pol_particle.theta_array(), 1e-8, 1e-5));
    assert!(tor_particle.zeta_array().relative_eq(&pol_particle.zeta_array(), 1e-8, 1e-5));
    assert!(tor_particle.rho_array().relative_eq(&pol_particle.rho_array(), 1e-8, 1e-5));
    assert!(tor_particle.ptheta_array().relative_eq(&pol_particle.ptheta_array(), 1e-8, 1e-5));
    assert!(tor_particle.pzeta_array().relative_eq(&pol_particle.pzeta_array(), 1e-8, 1e-5));
    assert!(tor_particle.energy_array().relative_eq(&pol_particle.energy_array(), 1e-8, 1e-5));
}
