import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version

LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)


def test_flute():
    mode = dex.FluteMode(1e-4, LCFS, 3, 2, 0)
    _test_mode_base(mode)
    assert isfinite(mode.psi_last)
    assert mode.lcfs.value == 0.05
    assert mode.phase == 0
    assert mode.epsilon == 1e-4


def test_nc(nc_flute_mode: dex.NcFluteMode):
    _test_mode_base(nc_flute_mode)
    assert isfinite(nc_flute_mode.psi_last)
    assert isfinite(nc_flute_mode.psip_last)
    assert nc_flute_mode.interp_type == "Cubic"
    assert isinstance(nc_flute_mode.path, str)
    assert isinstance(nc_flute_mode.netcdf_version, Version)
    assert nc_flute_mode.phase_method == "Interpolation"
    assert nc_flute_mode.analytical_threshold_index == 3
    assert nc_flute_mode.phase_average is None
    assert isinstance(nc_flute_mode.psi_array, np.ndarray)
    assert isinstance(nc_flute_mode.psip_array, np.ndarray)
    assert isinstance(nc_flute_mode.alpha_array, np.ndarray)
    assert isinstance(nc_flute_mode.phase_array, np.ndarray)


def _test_mode_base(mode: dex.ModeObject):

    mode.__repr__()
    mode.__str__()

    assert mode.machine_type in ["Numerical", "Analytical"]
    assert mode.psi_state in ["Good", "Bad"]
    assert mode.psip_state in ["Good", "Bad"]

    assert mode.m == 3
    assert mode.n == 2

    try:
        # 1 Parameter Evaluations
        methods = [
            mode.ampl_of_psi,
            mode.ampl_of_psip,
            mode.phase_of_psi,
            mode.phase_of_psip,
            mode.m_of_psi,
            mode.m_of_psip,
            mode.dm_dpsi,
            mode.dm_dpsip,
            mode.dm_of_psi_dtheta,
            mode.dm_of_psip_dtheta,
            mode.dm_of_psi_dzeta,
            mode.dm_of_psip_dzeta,
            mode.dm_of_psi_dt,
            mode.dm_of_psip_dt,
        ]

        # 0D evaluations
        flux = 1e-5
        theta = 1.57
        zeta = 1.57
        t = 0
        for method in methods:
            assert isfinite(method(flux, theta, zeta, t))
            assert isinstance(method(flux, theta, zeta, t), float)

        # 1D Evaluations
        fluxes = np.linspace(1e-5, 1e-4, 5)
        thetas = np.linspace(0, np.pi, 5)
        zetas = np.linspace(0, np.pi, 5)
        ts = np.linspace(0, 1, 5)
        for method in methods:
            assert method(fluxes, thetas, zetas, ts).ndim == 1
            assert isinstance(method(fluxes, thetas, zetas, ts), np.ndarray)

        # 4D Evaluations
        fluxes = np.random.random([2] * 4) * 1e-5
        thetas = np.random.random([2] * 4) * np.pi
        zetas = np.random.random([2] * 4) * np.pi
        ts = np.random.random([2] * 4) * 0
        assert fluxes.ndim == 4
        for method in methods:
            assert method(fluxes, thetas, zetas, ts).ndim == 4
            assert isinstance(method(fluxes, thetas, zetas, ts), np.ndarray)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )
