from pyncare import Qfactor, Currents, Bfield, Perturbation, Particle, MappingParameters


def test_particle_integration(
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
    perturbation: Perturbation,
    particle: Particle,
):
    assert particle.status == "Initialized"
    particle.integrate(
        qfactor=qfactor,
        currents=currents,
        bfield=bfield,
        perturbation=perturbation,
        t_eval=[0, 10],
    )
    assert particle.status == "Integrated"


def test_particle_mapping(
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
    perturbation: Perturbation,
    particle: Particle,
    params: MappingParameters,
):
    assert particle.status == "Initialized"
    particle.map(
        qfactor=qfactor,
        currents=currents,
        bfield=bfield,
        perturbation=perturbation,
        params=params,
    )
    assert particle.status == "Mapped"
