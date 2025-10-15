import poincare_maps as pm
import pytest


@pytest.fixture(scope="session")
def bfield():
    """Creates a Bfield object from "./data.nc" netCDF file."""
    return pm.Bfield("./data.nc", "bicubic")


@pytest.fixture(scope="session")
def qfactor():
    """Creates a Qfactor object from "./data.nc" netCDF file."""
    return pm.Qfactor("./data.nc", "akima")


@pytest.fixture(scope="session")
def current():
    """Creates a Current object from "./data.nc" netCDF file."""
    return pm.Current("./data.nc", "akima")


@pytest.fixture(scope="session")
def harmonic():
    """Creates a harmonic object."""
    return pm.Harmonic("./data.nc", "akima", 3, 2, 0)


@pytest.fixture(scope="session")
def perturbation():
    """Creates a perturbation object."""
    harmonics = [
        pm.Harmonic("./data.nc", "akima", 3, 2, 0),
        pm.Harmonic("./data.nc", "akima", 3, 1, 0),
    ]
    return pm.Perturbation(harmonics)


# ================================================================================================


@pytest.fixture(scope="session")
def initial():
    """Creates a set of initial conditions."""
    return pm.InitialConditions(t0=0, theta0=0, psip0=0.05, rho0=0.01, zeta0=0, mu=0)


@pytest.fixture(scope="function")
def particle(initial):
    """Creates a particle without running it."""
    return pm.Particle(initial)


@pytest.fixture(scope="function")
def poincare(particle):
    """Creates a poincare map object with some particles inside."""
    poincare = pm.Poincare()
    for _ in range(3):
        poincare.add_particle(particle)
    return poincare
