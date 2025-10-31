from pyncare import Qfactor, Current, Bfield, Perturbation, Particle


def test_particle_integration(
    qfactor: Qfactor,
    current: Current,
    bfield: Bfield,
    perturbation: Perturbation,
    particle: Particle,
):
    assert particle.status == "Initialized"
    particle.integrate(
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        per=perturbation,
        t_eval=[0, 10],
    )
    assert particle.status == "Integrated"
