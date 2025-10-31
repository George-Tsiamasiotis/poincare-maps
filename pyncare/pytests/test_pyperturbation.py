import pytest
from pyncare import Harmonic, Perturbation


def test_pyperturbation_derived_fields(perturbation: Perturbation):
    """Tests that the fields derived from the wrapped Rust object are the correct type."""
    assert isinstance(perturbation.harmonics, list)


def test_pypertyrbation_eval(perturbation: Perturbation):
    psip = 0.015
    theta = 1
    zeta = 2
    assert isinstance(perturbation.p(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dpsip(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dtheta(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dzeta(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dt(psip, theta, zeta), float)


def test_perturbation_getitem(perturbation: Perturbation):
    assert isinstance(perturbation[0], Harmonic)
    assert isinstance(perturbation[1], Harmonic)


def test_immutability(perturbation: Perturbation):
    """Tests that harmonic1 fields are immutable."""
    with pytest.raises(AttributeError):
        perturbation.harmonics = []


def test_repr(perturbation: Perturbation):
    str(perturbation)
