//! Tests `Equilibrium` functionality.

use dexter_equilibrium::*;

#[test]
#[rustfmt::skip]
fn create_analytical() {
    let geometry = LarGeometry::new(1.0, 1.75, 0.5);
    let psi_last = (geometry.rlast() / geometry.raxis()).powi(2) / 2.0;
    let lcfs = LastClosedFluxSurface::Toroidal(psi_last);

    let eq = dbg!(Equilibrium {
        geometry: Some(Box::new(LarGeometry::new(1.0, 1.75, 0.5))),
        qfactor: Box::new(ParabolicQfactor::new(1.1, 3.9, lcfs)),
        current: Box::new(LarCurrent::new()),
        bfield: Box::new(LarBfield::new()),
        perturbation: Perturbation::new(&[
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]),
    });

    let acc1 = &mut Accelerator::new();
    let acc2 = &mut Accelerator2d::new();
    let caches = &mut eq.perturbation.generate_caches();

    assert!(eq.psi_last().is_finite());
    assert!(eq.psip_last().is_finite());

    assert!(eq.geometry.is_some_and(|g| g.r_of_psi(0.01, acc1).unwrap().is_finite()));
    assert!(eq.qfactor.q_of_psi(0.01, acc1).unwrap().is_finite());
    assert!(eq.current.g_of_psi(0.01, acc1).unwrap().is_finite());
    assert!(eq.bfield.b_of_psi(0.01, 0.0, acc2).unwrap().is_finite());
    assert!(eq.perturbation.p_of_psi(0.01, 0.0, 0.0, 0.0, caches).unwrap().is_finite());
}
