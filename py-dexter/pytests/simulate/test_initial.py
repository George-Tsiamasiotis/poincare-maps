import pytest
import dexter as dex
from math import isclose


def test_initial_flux():
    with pytest.raises(RuntimeError):
        dex.InitialFlux()

    flux0 = dex.InitialFlux.Toroidal(0.01)
    assert flux0.kind == "Toroidal"
    assert isclose(flux0.value, 0.01)
    flux0.__repr__()
    flux0.__str__()

    flux0 = dex.InitialFlux.Poloidal(0.02)
    assert flux0.kind == "Poloidal"
    assert isclose(flux0.value, 0.02)
    flux0.__repr__()
    flux0.__str__()


def test_initial_conditions_boozer():
    with pytest.raises(RuntimeError):
        dex.InitialConditions()

    flux0 = dex.InitialFlux.Toroidal(0.01)
    init = dex.InitialConditions.Boozer(0, flux0, 1, 2, 3, 4)
    assert init.coordinate_set == "BoozerToroidal"
    assert isinstance(init.flux0, dex.InitialFlux)
    assert init.flux0.kind == "Toroidal"
    assert isclose(init.t0, 0)
    assert isclose(init.flux0.value, 0.01)
    assert isclose(init.theta0, 1)
    assert isclose(init.zeta0, 2)
    assert isclose(init.rho0, 3)
    assert isclose(init.mu0, 4)
    with pytest.raises(AttributeError):
        init.pzeta0
    init.__repr__()
    init.__str__()


def test_initial_conditions_mixed():
    flux0 = dex.InitialFlux.Poloidal(0.01)
    init = dex.InitialConditions.Mixed(0, flux0, 1, 2, 3, 4)
    assert init.coordinate_set == "MixedPoloidal"
    assert isinstance(init.flux0, dex.InitialFlux)
    assert init.flux0.kind == "Poloidal"
    assert isclose(init.t0, 0)
    assert isclose(init.flux0.value, 0.01)
    assert isclose(init.theta0, 1)
    assert isclose(init.zeta0, 2)
    assert isclose(init.pzeta0, 3)
    assert isclose(init.mu0, 4)
    with pytest.raises(AttributeError):
        init.rho0
    init.__repr__()
    init.__str__()
