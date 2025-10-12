import pytest
from pathlib import PosixPath


def test_bfield_fields(bfield):
    """Tests that bfield fields are the correct type."""
    assert isinstance(bfield.path, PosixPath)
    assert isinstance(bfield.typ, str)
    assert isinstance(bfield.baxis, float)
    assert isinstance(bfield.raxis, float)
    assert isinstance(bfield.psip_wall, float)
    assert isinstance(bfield.psi_wall, float)


def test_immutability(bfield):
    """Tests that bfield fields are immutable."""
    with pytest.raises(AttributeError):
        bfield.baxis += 1
        bfield.raxis += 1
        bfield.psip_wall += 1
        bfield.psi_wall += 1
        bfield.path = ""
        bfield.typ = ""


def test_data_extraction(bfield):
    """Tests that all extracted data are numpy arrays of the correct shape."""
    psip_data = bfield.psip_data()
    theta_data = bfield.theta_data()
    b_data = bfield.b_data()
    r_data = bfield.r_data()
    z_data = bfield.z_data()
    db_dpsip_data = bfield.db_dpsip_data()
    db_dtheta_data = bfield.db_dtheta_data()

    assert psip_data.shape == (101,)
    assert theta_data.shape == (3620,)
    assert b_data.shape == (101, 3620)
    assert r_data.shape == (101, 3620)
    assert z_data.shape == (101, 3620)
    assert db_dpsip_data.shape == (101, 3620)
    assert db_dtheta_data.shape == (101, 3620)


def test_bfield_print(bfield):
    """Tests __repr__() functionality."""
    str(bfield)
    print(bfield)
