import pytest
from pathlib import PosixPath
from pyncare import Qfactor


def test_pyqfactor_derived_fields(qfactor: Qfactor):
    """Tests that the fields derived from the wrapped Rust object are the correct type."""
    assert isinstance(qfactor.path, PosixPath)
    assert isinstance(qfactor.typ, str)
    assert isinstance(qfactor.psip_wall, float)
    assert isinstance(qfactor.psi_wall, float)


def test_pyqfactor_eval(qfactor: Qfactor):
    psip = 0.015
    assert isinstance(qfactor.q(psip), float)
    assert isinstance(qfactor.psi(psip), float)


def test_data_extraction(qfactor: Qfactor):
    """Tests that all extracted data are numpy arrays of the correct shape."""
    q_data = qfactor.q_data
    psip_data = qfactor.psip_data
    psi_data = qfactor.psi_data
    q_data_derived = qfactor.q_data_derived

    assert psip_data.ndim == 1
    assert q_data.ndim == 1
    assert psi_data.ndim == 1
    assert q_data_derived.ndim == 1


def test_immutability(qfactor: Qfactor):
    """Tests that qfactor fields are immutable."""
    with pytest.raises(AttributeError):
        qfactor.psip_wall += 1
        qfactor.psi_wall += 1
        qfactor.path = ""
        qfactor.typ = ""


def test_repr(qfactor: Qfactor):
    str(qfactor)
