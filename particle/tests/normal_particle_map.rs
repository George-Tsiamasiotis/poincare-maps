mod common;

use particle::*;

use crate::common::create_equilibrium;

#[test]
fn test_normal_particle_map_zeta() {
    let (qfactor, current, bfield, per) = create_equilibrium();
    let psip_wall = qfactor.psip_wall;

    let mapping = Mapping::new(PoincareSection::ConstZeta, 1.0, 6);
    let mut particle = Particle::new(0.0, 1.0, 0.5 * psip_wall, 0.01, 0.0, 0.0);
    assert!(matches!(particle.status, IntegrationStatus::Initialized));
    particle
        .map(&qfactor, &bfield, &current, &per, &mapping)
        .unwrap();
    dbg!(&particle);

    assert!(matches!(particle.status, IntegrationStatus::Integrated));
    assert_eq!(particle.evolution.zeta.len(), mapping.intersections + 1) // exclude initial point
}

#[test]
fn test_normal_particle_map_theta() {
    let (qfactor, current, bfield, per) = create_equilibrium();
    let psip_wall = qfactor.psip_wall;

    let mapping = Mapping::new(PoincareSection::ConstTheta, 1.0, 6);
    let mut particle = Particle::new(0.0, 0.0, 0.5 * psip_wall, 0.01, 1.0, 0.0);
    assert!(matches!(particle.status, IntegrationStatus::Initialized));
    particle
        .map(&qfactor, &bfield, &current, &per, &mapping)
        .unwrap();
    dbg!(&particle);

    assert!(matches!(particle.status, IntegrationStatus::Integrated));
    assert_eq!(particle.evolution.theta.len(), mapping.intersections + 1) // exclude initial point
}
