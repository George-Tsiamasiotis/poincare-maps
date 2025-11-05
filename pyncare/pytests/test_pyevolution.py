from pyncare._core import Particle


def test_evolution_fields(integrated_particle: Particle):
    assert isinstance(integrated_particle.evolution.steps_taken, int)
    assert isinstance(integrated_particle.evolution.steps_stored, int)


def test_evolution_extraction(integrated_particle: Particle):
    evolution = integrated_particle.evolution

    assert evolution.time.ndim == 1
    assert evolution.theta.ndim == 1
    assert evolution.psip.ndim == 1
    assert evolution.rho.ndim == 1
    assert evolution.zeta.ndim == 1
    assert evolution.psi.ndim == 1
    assert evolution.ptheta.ndim == 1
    assert evolution.pzeta.ndim == 1


def test_evolution_array_lengths(integrated_particle: Particle):
    evolution = integrated_particle.evolution

    length = len(evolution.time)

    assert len(evolution.theta) == length
    assert len(evolution.psip) == length
    assert len(evolution.rho) == length
    assert len(evolution.zeta) == length
    assert len(evolution.psi) == length
    assert len(evolution.ptheta) == length
    assert len(evolution.pzeta) == length


def test_repr(integrated_particle: Particle):
    str(integrated_particle.evolution)
