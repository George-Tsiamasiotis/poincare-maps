import pytest
from pathlib import PosixPath


def test_pybfield_derived_fields(bfield):
    """Tests that the fields derived from the wrapped Rust object are the correct type."""
    assert isinstance(bfield.path, PosixPath)
    assert isinstance(bfield.typ, str)
    assert isinstance(bfield.baxis, float)
    assert isinstance(bfield.raxis, float)
    assert isinstance(bfield.psip_wall, float)
    assert isinstance(bfield.psi_wall, float)


def test_pybfield_eval(bfield):
    psip = 0.015
    theta = 1
    assert isinstance(bfield.b(psip, theta), float)
    assert isinstance(bfield.db_dpsip(psip, theta), float)
    assert isinstance(bfield.db_dtheta(psip, theta), float)
    assert isinstance(bfield.d2b_dpsip2(psip, theta), float)
    assert isinstance(bfield.d2b_dtheta2(psip, theta), float)
    assert isinstance(bfield.d2b_dpsip_dtheta(psip, theta), float)


def test_data_extraction(bfield):
    """Tests that all extracted data are numpy arrays of the correct shape."""
    psip_data = bfield.psip_data
    theta_data = bfield.theta_data
    b_data = bfield.b_data
    r_data = bfield.r_data
    z_data = bfield.z_data
    db_dpsip_data = bfield.db_dpsip_data
    db_dtheta_data = bfield.db_dtheta_data

    assert psip_data.ndim == 1
    assert theta_data.ndim == 1
    assert b_data.ndim == 2
    assert r_data.ndim == 2
    assert z_data.ndim == 2
    assert db_dpsip_data.ndim == 2
    assert db_dtheta_data.ndim == 2


def test_immutability(bfield):
    """Tests that bfield fields are immutable."""
    with pytest.raises(AttributeError):
        bfield.psip_wall += 1
        bfield.psi_wall += 1
        bfield.baxis += 1
        bfield.raxis += 1
        bfield.path = ""
        bfield.typ = ""


def test_repr(bfield):
    str(bfield)
