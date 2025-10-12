import poincare_maps as pm
import pytest


@pytest.fixture(scope="session")
def bfield():
    """Creates a Bfield object from "./data.nc" netCDF file."""
    return pm.Bfield("./data.nc", "bicubic")
