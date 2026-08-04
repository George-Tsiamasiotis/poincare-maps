import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version

LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)


def test_unity():
    qfactor = dex.UnityQfactor(LCFS)
    assert qfactor.object_type == "Analytical"
    assert qfactor.psi_state == "Good"
    assert qfactor.psip_state == "Good"
    assert qfactor.psi_last == 0.05
    assert qfactor.psip_last == 0.05
    assert qfactor.qaxis == 1
    assert qfactor.qlast == 1
    qfactor.__repr__()
    qfactor.__str__()
    _test_qfactor_base(qfactor)
    with pytest.raises(Exception):
        qfactor.psi_of_q(1)
    with pytest.raises(Exception):
        qfactor.psip_of_q(1)


def test_parabolic():
    qfactor = dex.ParabolicQfactor(1.1, 3.9, LCFS)
    _test_qfactor_base(qfactor)
    assert qfactor.psi_last == 0.05
    assert qfactor.psip_last == qfactor.psip_of_psi(qfactor.psi_last)
    assert qfactor.qaxis == 1.1
    assert qfactor.qlast == 3.9


def test_nc(nc_qfactor: dex.NcQfactor):
    _test_qfactor_base(nc_qfactor)
    assert nc_qfactor.interp_type == "Cubic"
    assert isinstance(nc_qfactor.path, str)
    assert isinstance(nc_qfactor.netcdf_version, Version)
    assert isinstance(nc_qfactor.psi_array, np.ndarray)
    assert isinstance(nc_qfactor.psip_array, np.ndarray)
    assert isinstance(nc_qfactor.q_array, np.ndarray)


def _test_qfactor_base(qfactor: dex.QfactorObject):

    qfactor.__repr__()
    qfactor.__str__()

    assert qfactor.object_type in ["Numerical", "Analytical"]
    assert qfactor.psi_state in ["Good", "Bad"]
    assert qfactor.psip_state in ["Good", "Bad"]

    assert isfinite(qfactor.psi_last)
    assert isfinite(qfactor.psip_last)
    assert isfinite(qfactor.qlast)
    assert isfinite(qfactor.qaxis)

    methods = [
        qfactor.psip_of_psi,
        qfactor.psi_of_psip,
        qfactor.q_of_psi,
        qfactor.q_of_psip,
        qfactor.dpsip_dpsi,
        qfactor.dpsi_dpsip,
        qfactor.iota_of_psi,
        qfactor.iota_of_psip,
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
