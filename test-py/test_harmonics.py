import pytest
import numpy as np
from pathlib import PosixPath


def test_harmonic_fields(harmonic):
    assert isinstance(harmonic.path, PosixPath)
    assert isinstance(harmonic.typ, str)
    assert isinstance(harmonic.psip_wall, float)
    assert isinstance(harmonic.m, float)
    assert isinstance(harmonic.n, float)
    assert isinstance(harmonic.phase, float)
    assert isinstance(harmonic.amax, float)


def test_immutability(harmonic):
    """Tests that harmonic fields are immutable."""
    with pytest.raises(AttributeError):
        harmonic.amax += 1
        harmonic.phase += 1
        harmonic.m += 1
        harmonic.n += 1
        harmonic.psip_wall += 1
        harmonic.path = ""
        harmonic.typ = ""


def test_data_extraction(harmonic):
    """Tests that all extracted data are numpy arrays."""
    psip_data = harmonic.psip_data()
    a_data = harmonic.a_data()

    assert isinstance(psip_data, np.ndarray)
    assert isinstance(a_data, np.ndarray)


def test_harmonic_print(harmonic):
    """Tests __repr__() functionality."""
    str(harmonic)
    print(harmonic)
