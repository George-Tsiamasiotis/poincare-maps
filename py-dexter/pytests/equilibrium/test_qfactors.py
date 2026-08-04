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
    assert qfactor.object_type == "Analytical"
    assert qfactor.psi_state == "Good"
    assert qfactor.psip_state == "Good"
    assert qfactor.psi_last == 0.05
    assert qfactor.psip_last == qfactor.psip_of_psi(qfactor.psi_last)
    assert qfactor.qaxis == 1.1
    assert qfactor.qlast == 3.9
    qfactor.__repr__()
    qfactor.__str__()
    _test_qfactor_base(qfactor)


def test_nc(ncqfactor: dex.NcQfactor):
    assert ncqfactor.object_type == "Numerical"
    assert ncqfactor.psi_state == "Good"
    assert ncqfactor.psip_state == "Good"
    assert isfinite(ncqfactor.psi_last)
    assert isfinite(ncqfactor.psip_last)
    assert isfinite(ncqfactor.qlast)
    assert isfinite(ncqfactor.qaxis)
    ncqfactor.__repr__()
    ncqfactor.__str__()
    _test_qfactor_base(ncqfactor)
    assert ncqfactor.interp_type == "Cubic"
    assert isinstance(ncqfactor.path, str)
    assert isinstance(ncqfactor.netcdf_version, Version)
    assert isinstance(ncqfactor.psi_array, np.ndarray)
    assert isinstance(ncqfactor.psip_array, np.ndarray)
    assert isinstance(ncqfactor.q_array, np.ndarray)


def _test_qfactor_base(qfactor: dex.QfactorObject):

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
