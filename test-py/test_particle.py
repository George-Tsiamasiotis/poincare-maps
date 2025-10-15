def test_particle_run_ode(bfield, qfactor, current, perturbation, particle):
    str(particle)
    particle.run_ode(
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        per=perturbation,
        t_eval=(0, 20),
        steps=0,
    )
    str(particle)


def test_particle_run_henon_theta(bfield, qfactor, current, perturbation, particle):
    str(particle)
    particle.run_henon(
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        per=perturbation,
        intersection=3.14,
        angle="theta",
        turns=3,
    )
    str(particle)


def test_particle_run_henon_zeta(bfield, qfactor, current, perturbation, particle):
    str(particle)
    particle.run_henon(
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        per=perturbation,
        intersection=3.14,
        angle="zeta",
        turns=3,
    )
    str(particle)
