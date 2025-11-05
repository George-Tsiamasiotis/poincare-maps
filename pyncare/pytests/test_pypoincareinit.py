from pyncare import PoincareInit
import numpy as np
import pytest


def test_poincare_init_fields(poincare_init: PoincareInit):
    assert isinstance(poincare_init.thetas, np.ndarray)
    assert isinstance(poincare_init.psips, np.ndarray)
    assert isinstance(poincare_init.rhos, np.ndarray)
    assert isinstance(poincare_init.zetas, np.ndarray)
    assert isinstance(poincare_init.mus, np.ndarray)


def test_immutability(poincare_init: PoincareInit):
    with pytest.raises(AttributeError):
        poincare_init.thetas += 1
        poincare_init.psips += 1
        poincare_init.rhos += 1
        poincare_init.zetas += 1
        poincare_init.mus += 1


def test_repr(poincare_init: PoincareInit):
    str(poincare_init)
