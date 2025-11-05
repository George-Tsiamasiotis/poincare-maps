import pyncare as pc
import pytest
import numpy as np


@pytest.fixture(scope="session")
def qfactor():
    """Creates a Qfactor object from "./data.nc" netCDF file."""
    return pc.Qfactor("./data.nc", "akima")


@pytest.fixture(scope="session")
def current():
    """Creates a Current object from "./data.nc" netCDF file."""
    return pc.Current("./data.nc", "akima")


@pytest.fixture(scope="session")
def bfield():
    """Creates a Bfield object from "./data.nc" netCDF file."""
    return pc.Bfield("./data.nc", "bicubic")


@pytest.fixture(scope="session")
def harmonic1():
    """Creates a Harmonic object from "./data.nc" netCDF file."""
    return pc.Harmonic("./data.nc", "akima", m=1, n=2, phase=0)


@pytest.fixture(scope="session")
def harmonic2():
    """Creates a Harmonic object from "./data.nc" netCDF file."""
    return pc.Harmonic("./data.nc", "akima", m=3, n=2, phase=0)


@pytest.fixture(scope="session")
def perturbation(harmonic1, harmonic2):
    """Creates a Perturbation object with the 2 fixture harmonics."""
    return pc.Perturbation(harmonics=[harmonic1, harmonic2])


# =========================================================================================


@pytest.fixture(scope="session")
def initial_conditions(qfactor):
    """Creates an InitialConditons object."""
    psip_wall = qfactor.psip_wall
    return pc.InitialConditions(
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
    return pc.Particle(initial_conditions)


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
    return pc.Mapping("theta", 3.14, 10)


# =========================================================================================


@pytest.fixture(scope="session")
def poincare_init(qfactor: pc.Qfactor):
    """Creates a PoincareInit object."""
    num = 5
    return pc.PoincareInit(
        thetas=np.linspace(0, np.pi, num),
        psips=np.linspace(0.1, 0.7, num) * qfactor.psip_wall,
        rhos=0.001 * np.ones(num),
        zetas=np.zeros(num),
        mus=np.zeros(num),
    )


@pytest.fixture(scope="session")
def poincare(poincare_init: pc.PoincareInit, mapping: pc.Mapping):
    """Creates a Poincare object."""
    return pc.Poincare(init=poincare_init, mapping=mapping)
