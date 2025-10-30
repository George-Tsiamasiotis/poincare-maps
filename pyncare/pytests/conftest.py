import pyncare as pm
import pytest


@pytest.fixture(scope="session")
def qfactor():
    """Creates a Qfactor object from "./data.nc" netCDF file."""
    return pm.Qfactor("./data.nc", "akima")


@pytest.fixture(scope="session")
def current():
    """Creates a Current object from "./data.nc" netCDF file."""
    return pm.Current("./data.nc", "akima")


@pytest.fixture(scope="session")
def bfield():
    """Creates a Bfield object from "./data.nc" netCDF file."""
    return pm.Bfield("./data.nc", "bicubic")


@pytest.fixture(scope="session")
def harmonic1():
    """Creates a Harmonic object from "./data.nc" netCDF file."""
    return pm.Harmonic("./data.nc", "akima", m=1, n=2, phase=0)


@pytest.fixture(scope="session")
def harmonic2():
    """Creates a Harmonic object from "./data.nc" netCDF file."""
    return pm.Harmonic("./data.nc", "akima", m=3, n=2, phase=0)
