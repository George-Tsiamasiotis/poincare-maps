import pytest
from pyncare import Harmonic


def test_pyharmonic1_derived_fields(harmonic1: Harmonic):
    assert isinstance(harmonic1.path, str)
    assert isinstance(harmonic1.typ, str)
    assert isinstance(harmonic1.m, int)
    assert isinstance(harmonic1.n, int)
    assert isinstance(harmonic1.phase_average, float)
    assert isinstance(harmonic1.psip_wall, float)


def test_pyharmonic1_eval(harmonic1: Harmonic):
    psip = 0.015
    theta = 1
    zeta = 2
    assert isinstance(harmonic1.h(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dpsip(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dtheta(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dzeta(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dt(psip, theta, zeta), float)


def test_data_extraction(harmonic1: Harmonic):
    assert harmonic1.psip_data.ndim == 1
    assert harmonic1.a_data.ndim == 1


def test_immutability(harmonic1: Harmonic):
    with pytest.raises(AttributeError):
        harmonic1.psip_wall += 1
        harmonic1.m += 1
        harmonic1.n += 1
        harmonic1.phase_average += 1
        harmonic1.path = ""
        harmonic1.typ = ""


def test_repr(harmonic1: Harmonic):
    str(harmonic1)
