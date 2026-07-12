//! Common testing utilities.

use dexter_equilibrium::*;
use dexter_simulate::*;

#[allow(dead_code, reason = "used in tests")]
pub(crate) fn lar_equilibrium() -> Equilibrium {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    Equilibrium {
        geometry: None,
        qfactor: Box::new(ParabolicQfactor::new(1.1, 3.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(&[
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 3, 1, 0.0)),
        ]),
    }
}

/// Checks that the time series calculated from an integration routine are valid.
#[allow(dead_code, reason = "used in tests")]
pub(crate) fn check_integrated_particle_arrays(particle: &Particle) {
    let _: InitialConditions = particle.initial_conditions();
    assert_eq!(particle.t_array().len(), particle.steps_stored());
    assert_eq!(particle.psi_array().len(), particle.steps_stored());
    assert_eq!(particle.psip_array().len(), particle.steps_stored());
    assert_eq!(particle.theta_array().len(), particle.steps_stored());
    assert_eq!(particle.zeta_array().len(), particle.steps_stored());
    assert_eq!(particle.rho_array().len(), particle.steps_stored());
    assert_eq!(particle.mu_array().len(), particle.steps_stored());
    assert_eq!(particle.ptheta_array().len(), particle.steps_stored());
    assert_eq!(particle.pzeta_array().len(), particle.steps_stored());
    assert_eq!(particle.energy_array().len(), particle.steps_stored());

    assert!(particle.t_array().iter().all(|x| x.is_finite()));
    assert!(particle.psi_array().iter().all(|x| x.is_finite()));
    assert!(particle.psip_array().iter().all(|x| x.is_finite()));
    assert!(particle.theta_array().iter().all(|x| x.is_finite()));
    assert!(particle.zeta_array().iter().all(|x| x.is_finite()));
    assert!(particle.rho_array().iter().all(|x| x.is_finite()));
    assert!(particle.mu_array().iter().all(|x| x.is_finite()));
    assert!(particle.ptheta_array().iter().all(|x| x.is_finite()));
    assert!(particle.pzeta_array().iter().all(|x| x.is_finite()));
    assert!(particle.energy_array().iter().all(|x| x.is_finite()));
}
