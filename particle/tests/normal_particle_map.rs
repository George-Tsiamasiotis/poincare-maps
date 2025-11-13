mod common;

use particle::*;

use crate::common::create_equilibrium;

#[test]
fn test_normal_particle_map_zeta() {
    let (qfactor, currents, bfield, perturbation) = create_equilibrium();
    let psip_wall = qfactor.psip_wall();

    let intial = InitialConditions {
        time0: 0.0,
        theta0: 1.0,
        psip0: 0.5 * psip_wall,
        rho0: 0.0001,
        zeta0: 0.0,
        mu: 0.0,
    };

    let params = MappingParameters::new(PoincareSection::ConstZeta, 1.0, 6);
    let mut particle = Particle::new(&intial);
    assert!(matches!(particle.status, IntegrationStatus::Initialized));
    particle
        .map(&qfactor, &bfield, &currents, &perturbation, &params)
        .unwrap();
    dbg!(&particle);

    assert!(matches!(particle.status, IntegrationStatus::Mapped));
    assert_eq!(particle.evolution.zeta.len(), params.intersections + 1) // exclude initial point
}

#[test]
fn test_normal_particle_map_theta() {
    let (qfactor, currents, bfield, perturbation) = create_equilibrium();
    let psip_wall = qfactor.psip_wall();

    let intial = InitialConditions {
        time0: 0.0,
        theta0: 1.0,
        psip0: 0.5 * psip_wall,
        rho0: 0.0001,
        zeta0: 0.0,
        mu: 0.0,
    };

    let params = MappingParameters::new(PoincareSection::ConstTheta, 1.0, 6);
    let mut particle = Particle::new(&intial);
    assert!(matches!(particle.status, IntegrationStatus::Initialized));
    particle
        .map(&qfactor, &bfield, &currents, &perturbation, &params)
        .unwrap();
    dbg!(&particle);

    assert!(matches!(particle.status, IntegrationStatus::Mapped));
    assert_eq!(particle.evolution.theta.len(), params.intersections + 1) // exclude initial point
}
