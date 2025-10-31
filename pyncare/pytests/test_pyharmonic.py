import pytest
from pathlib import PosixPath
from pyncare import Harmonic


def test_pyharmonic1_derived_fields(harmonic1: Harmonic):
    """Tests that the fields derived from the wrapped Rust object are the correct type."""
    assert isinstance(harmonic1.path, PosixPath)
    assert isinstance(harmonic1.typ, str)
    assert isinstance(harmonic1.m, float)
    assert isinstance(harmonic1.n, float)
    assert isinstance(harmonic1.phase, float)
    assert isinstance(harmonic1.amax, float)
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
    """Tests that all extracted data are numpy arrays of the correct shape."""
    psip_data = harmonic1.psip_data
    a_data = harmonic1.a_data

    assert psip_data.ndim == 1
    assert a_data.ndim == 1


def test_immutability(harmonic1: Harmonic):
    """Tests that harmonic1 fields are immutable."""
    with pytest.raises(AttributeError):
        harmonic1.psip_wall += 1
        harmonic1.m += 1
        harmonic1.n += 1
        harmonic1.phase += 1
        harmonic1.amax += 1
        harmonic1.path = ""
        harmonic1.typ = ""


def test_repr(harmonic1: Harmonic):
    str(harmonic1)
