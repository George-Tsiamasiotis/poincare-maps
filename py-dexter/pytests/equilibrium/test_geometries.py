import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version


def test_lar():
    geometry = dex.LarGeometry(2, 1.75, 0.5)
    _test_geometry_base(geometry)


def test_nc(nc_geometry: dex.NcGeometry):
    nc_geometry.__repr__()
    nc_geometry.__str__()
    _test_geometry_base(nc_geometry)
    assert isinstance(nc_geometry.path, str)
    assert isinstance(nc_geometry.netcdf_version, Version)
    assert nc_geometry.interp1d_type == "Cubic"
    assert nc_geometry.interp2d_type == "Bicubic"
    assert isfinite(nc_geometry.psi_last)
    assert isfinite(nc_geometry.psip_last)
    assert isinstance(nc_geometry.psi_array, np.ndarray)
    assert isinstance(nc_geometry.psip_array, np.ndarray)
    assert isinstance(nc_geometry.theta_array, np.ndarray)
    assert isinstance(nc_geometry.r_array, np.ndarray)
    rlab_array = nc_geometry.rlab_array
    assert isinstance(rlab_array, np.ndarray)
    assert rlab_array.shape == nc_geometry.shape
    zlab_array = nc_geometry.zlab_array
    assert isinstance(zlab_array, np.ndarray)
    assert zlab_array.shape == nc_geometry.shape
    jacobian_array = nc_geometry.jacobian_array
    assert isinstance(jacobian_array, np.ndarray)
    assert jacobian_array.shape == nc_geometry.shape


def _test_geometry_base(geometry: dex.GeometryObject):

    geometry.__repr__()
    geometry.__str__()

    assert geometry.object_type in ["Numerical", "Analytical"]
    assert geometry.psi_state in ["Good", "Bad"]
    assert geometry.psip_state in ["Good", "Bad"]

    assert isfinite(geometry.baxis)
    assert isfinite(geometry.raxis)
    assert isfinite(geometry.zaxis)
    assert isfinite(geometry.rgeo)
    assert isfinite(geometry.rlast)
    assert isinstance(geometry.rlab_last, np.ndarray)
    assert isinstance(geometry.zlab_last, np.ndarray)

    # 1 Parameter Evaluations
    methods = [
        geometry.r_of_psi,
        geometry.r_of_psip,
        geometry.psi_of_r,
        geometry.psip_of_r,
    ]

    try:
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
        fluxes = np.random.random([2] * 4) * 1e-5
        assert fluxes.ndim == 4
        for method in methods:
            assert method(fluxes).ndim == 4
            assert isinstance(method(fluxes), np.ndarray)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )

    # 2 Parameter Evaluations
    methods = [
        geometry.rlab_of_psi,
        geometry.rlab_of_psip,
        geometry.zlab_of_psi,
        geometry.zlab_of_psip,
        geometry.jacobian_of_psi,
        geometry.jacobian_of_psip,
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
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )
