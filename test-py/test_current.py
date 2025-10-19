import pytest
from pathlib import PosixPath


def test_current_fields(current):
    """Tests that current fields are the correct type."""
    assert isinstance(current.path, PosixPath)
    assert isinstance(current.typ, str)
    assert isinstance(current.psip_wall, float)
    assert isinstance(current.psi_wall, float)


def test_immutability(current):
    """Tests that current fields are immutable."""
    with pytest.raises(AttributeError):
        current.psip_wall += 1
        current.psi_wall += 1
        current.path = ""
        current.typ = ""


def test_data_extraction(current):
    """Tests that all extracted data are numpy arrays of the correct shape."""
    psip_data = current.psip_data()
    g_data = current.g_data()
    i_data = current.i_data()
    dg_dpsip_data = current.dg_dpsip_data()
    di_dpsip_data = current.di_dpsip_data()

    assert psip_data.ndim == 1
    assert g_data.ndim == 1
    assert i_data.ndim == 1
    assert dg_dpsip_data.ndim == 1
    assert di_dpsip_data.ndim == 1


def test_current_print(current):
    """Tests __repr__() functionality."""
    str(current)
    print(current)
