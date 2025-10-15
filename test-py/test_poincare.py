def test_poincare_run_theta(bfield, qfactor, current, perturbation, poincare):
    poincare.get_particles()
    poincare.get_angles()
    poincare.get_fluxes()
    poincare.run(
        qfactor=qfactor,
        bfield=bfield,
        current=current,
        per=perturbation,
        angle="zeta",
        intersection=3.14,
        turns=3,
    )
    poincare.get_particles()
    poincare.get_angles()
    poincare.get_fluxes()
