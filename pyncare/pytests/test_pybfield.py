import pytest
from pyncare import Bfield


def test_pybfield_derived_fields(bfield: Bfield):
    assert isinstance(bfield.path, str)
    assert isinstance(bfield.typ, str)
    assert isinstance(bfield.baxis, float)
    # assert isinstance(bfield.raxis, float)
    assert isinstance(bfield.psip_wall, float)


def test_pybfield_eval(bfield: Bfield):
    psip = 0.015
    theta = 1
    assert isinstance(bfield.b(psip, theta), float)
    assert isinstance(bfield.db_dpsip(psip, theta), float)
    assert isinstance(bfield.db_dtheta(psip, theta), float)
    assert isinstance(bfield.d2b_dpsip2(psip, theta), float)
    assert isinstance(bfield.d2b_dtheta2(psip, theta), float)
    assert isinstance(bfield.d2b_dpsip_dtheta(psip, theta), float)


def test_data_extraction(bfield: Bfield):
    assert bfield.psip_data.ndim == 1
    assert bfield.theta_data.ndim == 1
    assert bfield.b_data.ndim == 2
    assert bfield.r_data.ndim == 2
    assert bfield.z_data.ndim == 2
    assert bfield.db_dpsip_data.ndim == 2
    assert bfield.db_dtheta_data.ndim == 2


def test_immutability(bfield: Bfield):
    with pytest.raises(AttributeError):
        bfield.psip_wall += 1
        bfield.baxis += 1
        bfield.raxis += 1
        bfield.path = ""
        bfield.typ = ""


def test_repr(bfield: Bfield):
    str(bfield)
