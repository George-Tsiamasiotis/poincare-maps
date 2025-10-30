mod common;

use particle::*;

use crate::common::create_equilibrium;

#[test]
fn test_normal_particle_int() {
    let (qfactor, current, bfield, per) = create_equilibrium();
    let psip_wall = qfactor.psip_wall;

    let initial = InitialConditions {
        t0: 0.0,
        theta0: 1.0,
        psip0: 0.5 * psip_wall,
        rho0: 0.0001,
        zeta0: 0.0,
        mu: 0.0,
    };

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.status, IntegrationStatus::Initialized));
    particle
        .integrate(&qfactor, &bfield, &current, &per, (0.0, 500000.0))
        .unwrap();

    let hcache = particle.final_state.hcache.first().unwrap();

    assert!(matches!(particle.status, IntegrationStatus::Integrated));
    assert!(hcache.hits as f64 / hcache.misses as f64 - 3.0 <= 1e-5);

    dbg!(&particle);
}
