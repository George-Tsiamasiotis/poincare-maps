import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version


def test_lar():
    bfield = dex.LarBfield()
    _test_bfield_base(bfield)


def test_nc(nc_bfield: dex.NcBfield):
    _test_bfield_base(nc_bfield)
    assert isinstance(nc_bfield.path, str)
    assert isinstance(nc_bfield.netcdf_version, Version)
    assert nc_bfield.interp_type == "Bicubic"
    assert isfinite(nc_bfield.baxis)
    assert isfinite(nc_bfield.padding)
    assert isfinite(nc_bfield.psi_last)
    assert isfinite(nc_bfield.psip_last)
    assert isinstance(nc_bfield.psi_array, np.ndarray)
    assert isinstance(nc_bfield.psip_array, np.ndarray)
    assert isinstance(nc_bfield.theta_array, np.ndarray)
    assert isinstance(nc_bfield.theta_array_padded, np.ndarray)
    b_array = nc_bfield.b_array
    assert isinstance(b_array, np.ndarray)
    assert b_array.shape == nc_bfield.shape
    b_array_padded = nc_bfield.b_array_padded
    assert isinstance(b_array_padded, np.ndarray)
    assert b_array_padded.shape == nc_bfield.shape_padded


def _test_bfield_base(bfield: dex.BfieldObject):

    bfield.__repr__()
    bfield.__str__()

    assert bfield.machine_type in ["Numerical", "Analytical"]
    assert bfield.psi_state in ["Good", "Bad"]
    assert bfield.psip_state in ["Good", "Bad"]

    methods = [
        bfield.b_of_psi,
        bfield.b_of_psip,
        bfield.db_dpsi,
        bfield.db_dpsip,
        bfield.db_of_psi_dtheta,
        bfield.db_of_psip_dtheta,
    ]

    try:

        # 0D evaluations
        flux = 1e-5
        theta = 1.57
        for method in methods:
            assert isfinite(method(flux, theta))
            assert isinstance(method(flux, theta), float)

        # 1D Evaluations
        fluxes = np.linspace(1e-5, 1e-4, 5)
        thetas = np.linspace(0, np.pi, 5)
        for method in methods:
            assert method(fluxes, thetas).ndim == 1
            assert isinstance(method(fluxes, thetas), np.ndarray)

        # 4D Evaluations
        fluxes = np.random.random([2] * 4) * 1e-5
        thetas = np.random.random([2] * 4) * np.pi
        assert fluxes.ndim == 4
        for method in methods:
            assert method(fluxes, thetas).ndim == 4
            assert isinstance(method(fluxes, thetas), np.ndarray)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError("only testing the vectorized functions here")
