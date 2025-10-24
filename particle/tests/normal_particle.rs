mod common;

use particle::*;

use crate::common::create_equilibrium;

#[test]
fn test_normal_particle() {
    let (qfactor, current, bfield, per) = create_equilibrium();
    let psip_wall = qfactor.psip_wall;

    let mut particle = Particle::new(0.0, 1.0, 0.5 * psip_wall, 0.0001, 0.0, 0.0);
    assert!(matches!(particle.status, IntegrationStatus::Initialized));
    particle
        .integrate(&qfactor, &bfield, &current, &per, (0.0, 500000.0))
        .unwrap();

    assert!(matches!(particle.status, IntegrationStatus::Integrated));
    dbg!(&particle);
}
