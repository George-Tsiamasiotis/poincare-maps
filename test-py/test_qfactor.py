import pytest
from pathlib import PosixPath


def test_qfactor_fields(qfactor):
    """Tests that qfactor fields are the correct type."""
    assert isinstance(qfactor.path, PosixPath)
    assert isinstance(qfactor.typ, str)
    assert isinstance(qfactor.psip_wall, float)
    assert isinstance(qfactor.psi_wall, float)


def test_immutability(qfactor):
    """Tests that qfactor fields are immutable."""
    with pytest.raises(AttributeError):
        qfactor.psip_wall += 1
        qfactor.psi_wall += 1
        qfactor.path = ""
        qfactor.typ = ""


def test_data_extraction(qfactor):
    """Tests that all extracted data are numpy arrays of the correct shape."""
    psip_data = qfactor.psip_data()
    q_data = qfactor.psip_data()
    psi_data = qfactor.psip_data()

    assert psip_data.ndim == 1
    assert q_data.ndim == 1
    assert psi_data.ndim == 1


def test_qfactor_print(qfactor):
    """Tests __repr__() functionality."""
    str(qfactor)
    print(qfactor)
