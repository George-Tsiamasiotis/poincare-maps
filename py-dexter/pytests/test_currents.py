import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version


def test_lar():
    current = dex.LarCurrent()
    assert current.object_type == "Analytical"
    assert current.psi_state == "Good"
    assert current.psip_state == "Good"
    current.__repr__()
    current.__str__()
    _test_current_base(current)


def test_nc(nc_current: dex.NcCurrent):
    assert nc_current.object_type == "Numerical"
    assert nc_current.psi_state == "Good"
    assert nc_current.psip_state == "Good"
    assert isfinite(nc_current.psi_last)
    assert isfinite(nc_current.psip_last)
    nc_current.__repr__()
    nc_current.__str__()
    _test_current_base(nc_current)
    assert nc_current.interp_type == "Cubic"
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.netcdf_version, Version)
    assert isinstance(nc_current.psi_array, np.ndarray)
    assert isinstance(nc_current.psip_array, np.ndarray)
    assert isinstance(nc_current.i_array, np.ndarray)
    assert isinstance(nc_current.g_array, np.ndarray)


def _test_current_base(current: dex.CurrentObject):

    methods = [
        current.g_of_psi,
        current.g_of_psip,
        current.i_of_psi,
        current.i_of_psip,
        current.dg_dpsi,
        current.dg_dpsip,
        current.di_dpsi,
        current.di_dpsip,
    ]

    # 0D evaluations
    flux = 1e-5
    for method in methods:
        assert isfinite(method(flux))
        assert isinstance(method(flux), float)

    # 1D Evaluations
    fluxes = np.linspace(1e-5, 1e-4, 5)
    for method in methods:
        assert method(fluxes).ndim == 1
        assert isinstance(method(fluxes), np.ndarray)

    # 4D Evaluations
    grid = np.random.random([2] * 4) * 1e-5
    assert grid.ndim == 4
    for method in methods:
        assert method(grid).ndim == 4
        assert isinstance(method(grid), np.ndarray)
