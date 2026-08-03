import dexter as dex
import pytest


def test_lcfs():
    with pytest.raises(RuntimeError):
        dex.LastClosedFluxSurface()

    lcfs = dex.LastClosedFluxSurface.Toroidal(0.1)
    assert lcfs.value == 0.1
    assert lcfs.kind == "Toroidal"
    lcfs = dex.LastClosedFluxSurface.Poloidal(0.2)
    assert lcfs.value == 0.2
    assert lcfs.kind == "Poloidal"
