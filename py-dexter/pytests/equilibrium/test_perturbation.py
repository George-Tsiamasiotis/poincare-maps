import numpy as np
import dexter as dex

from math import isfinite

LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)


def test_empty_perturbation():
    per = dex.Perturbation()
    assert len(per) == 0
    per.__repr__()
    per.__str__()
    _test_evals(per)


def test_flute_perturbation():
    mode1 = dex.FluteMode(1e-4, LCFS, 3, 2, 0)
    mode2 = dex.FluteMode(2e-4, LCFS, 4, 3, 0)
    mode3 = dex.FluteMode(3e-4, LCFS, 5, 4, 0)
    per = dex.Perturbation([mode1, mode2, mode3])
    assert len(per) == 3
    per.__repr__()
    per.__str__()
    _test_evals(per)


def test_nc_perturbation(nc_flute_mode: dex.NcFluteMode):
    per = dex.Perturbation([nc_flute_mode, nc_flute_mode])
    assert len(per) == 2
    per.__repr__()
    per.__str__()
    _test_evals(per)


def _test_evals(per: dex.Perturbation):
    try:
        methods = [
            per.p_of_psi,
            per.p_of_psip,
            per.dp_dpsi,
            per.dp_dpsip,
            per.dp_of_psi_dtheta,
            per.dp_of_psip_dtheta,
            per.dp_of_psi_dzeta,
            per.dp_of_psip_dzeta,
            per.dp_of_psi_dt,
            per.dp_of_psip_dt,
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
