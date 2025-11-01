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


@pytest.fixture(scope="session")
def perturbation(harmonic1, harmonic2):
    """Creates a Perturbation object with the 2 fixture harmonics."""
    return pm.Perturbation(harmonics=[harmonic1, harmonic2])


@pytest.fixture(scope="session")
def initial_conditions(qfactor):
    """Creates an InitialConditons object."""
    psip_wall = qfactor.psip_wall
    return pm.InitialConditions(
        t0=0,
        theta0=3.14,
        psip0=0.5 * psip_wall,
        rho0=0.001,
        zeta0=0.0,
        mu=0,
    )


@pytest.fixture(scope="function")
def particle(initial_conditions):
    """Creates a Particle object from the initial conditions fixture."""
    return pm.Particle(initial_conditions)


@pytest.fixture(scope="function")
def integrated_particle(qfactor, current, bfield, perturbation, particle):
    """Creates a Particle object and integrates it."""
    particle.integrate(
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        per=perturbation,
        t_eval=[0, 10],
    )
    return particle


@pytest.fixture(scope="session")
def mapping():
    """Creates a Mapping object."""
    return pm.Mapping("theta", 3.14, 10)
