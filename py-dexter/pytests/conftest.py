import numpy
import pytest
import matplotlib
import dexter as dex

from math import sqrt

TEST_NETCDF_PATH = "./crates/dexter-machine/test_netcdf.nc"
TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-machine/toroidal_test_netcdf.nc"
POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-machine/poloidal_test_netcdf.nc"


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


@pytest.fixture(scope="session")
def nc_geometry() -> dex.NcGeometry:
    """Creates an NcGeometry object from the test NetCDF file."""
    return dex.NcGeometry(TEST_NETCDF_PATH, "Cubic", "Bicubic")


@pytest.fixture(scope="session")
def nc_flute_mode() -> dex.NcFluteMode:
    """Creates an NcFluteMode object from the test NetCDF file."""
    return dex.NcFluteMode(TEST_NETCDF_PATH, "Cubic", 3, 2)


@pytest.fixture(scope="session")
def lar_machine() -> dex.Machine:
    """Creates a typical LAR configuration."""
    LCFS = dex.LastClosedFluxSurface.Toroidal(0.03)
    raxis = 1.75
    rlast = sqrt(2 * LCFS.value) * raxis  # [m]
    return dex.Machine(
        geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
        qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
        current=dex.LarCurrent(),
        bfield=dex.LarBfield(),
    )


@pytest.fixture(scope="session")
def lar_machine_perturbed(lar_machine: dex.Machine) -> dex.Machine:
    """Creates a typical LAR configuration with two flute modes."""
    LCFS = dex.LastClosedFluxSurface.Toroidal(0.03)
    raxis = 1.75
    rlast = sqrt(2 * LCFS.value) * raxis  # [m]
    return dex.Machine(
        geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
        qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
        current=dex.LarCurrent(),
        bfield=dex.LarBfield(),
        perturbation=dex.Perturbation(
            [
                dex.FluteMode(1e-4, LCFS, 3, 1, 0),
                dex.FluteMode(1e-4, LCFS, 3, 2, 0),
            ]
        ),
    )
