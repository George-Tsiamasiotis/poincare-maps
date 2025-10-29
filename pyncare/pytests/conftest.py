import pyncare as pm
import pytest


@pytest.fixture(scope="session")
def qfactor():
    """Creates a Qfactor object from "./data.nc" netCDF file."""
    return pm.Qfactor("./data.nc", "akima")
