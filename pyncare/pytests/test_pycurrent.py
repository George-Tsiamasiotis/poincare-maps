import pytest
from pyncare import Current


def test_pycurrent_derived_fields(current: Current):
    assert isinstance(current.path, str)
    assert isinstance(current.typ, str)
    assert isinstance(current.psip_wall, float)


def test_pycurrent_eval(current: Current):
    psip = 0.015
    assert isinstance(current.g(psip), float)
    assert isinstance(current.i(psip), float)
    assert isinstance(current.dg_dpsip(psip), float)
    assert isinstance(current.di_dpsip(psip), float)


def test_data_extraction(current: Current):
    assert current.psip_data.ndim == 1
    assert current.g_data.ndim == 1
    assert current.i_data.ndim == 1


def test_immutability(current: Current):
    """Tests that current fields are immutable."""
    with pytest.raises(AttributeError):
        current.psip_wall += 1
        current.path = ""
        current.typ = ""


def test_repr(current):
    str(current)
