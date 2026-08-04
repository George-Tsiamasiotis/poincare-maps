import numpy
import pytest
import matplotlib
import dexter as dex

TEST_NETCDF_PATH = "./crates/dexter-equilibrium/test_netcdf.nc"
TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/toroidal_test_netcdf.nc"
POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/poloidal_test_netcdf.nc"


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace):
    matplotlib.use("agg")  # Disable interactive plots
    doctest_namespace["path"] = TEST_NETCDF_PATH
    doctest_namespace["dex"] = dex
    doctest_namespace["np"] = numpy


@pytest.fixture(scope="session")
def nc_qfactor() -> dex.NcQfactor:
    """Creates an NcQfactor object from the test NetCDF file."""
    return dex.NcQfactor(TEST_NETCDF_PATH, "Cubic")


@pytest.fixture(scope="session")
def nc_current() -> dex.NcCurrent:
    """Creates an NcCurrent object from the test NetCDF file."""
    return dex.NcCurrent(TEST_NETCDF_PATH, "Cubic")


@pytest.fixture(scope="session")
def nc_bfield() -> dex.NcBfield:
    """Creates an NcBfield object from the test NetCDF file."""
    return dex.NcBfield(TEST_NETCDF_PATH, "Bicubic")
