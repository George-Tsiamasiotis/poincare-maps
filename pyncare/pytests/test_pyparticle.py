from pyncare import Qfactor, Current, Bfield, Perturbation, Particle, Mapping


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


def test_particle_mapping(
    qfactor: Qfactor,
    current: Current,
    bfield: Bfield,
    perturbation: Perturbation,
    particle: Particle,
    mapping: Mapping,
):
    assert particle.status == "Initialized"
    particle.map(
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        per=perturbation,
        mapping=mapping,
    )
    assert particle.status == "Integrated"
