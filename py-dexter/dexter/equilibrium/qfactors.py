"""Defines the different q-factor objects as wrappers over `_PyQfactor`."""

import numpy as np

from dexter._core import _PyQfactor

from .utils import LastClosedFluxSurface
from .base import EquilibriumObject, FluxCommute, Qfactor


class UnityQfactor(EquilibriumObject, FluxCommute, Qfactor):
    r"""Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.

    Parameters
    ----------
    lcfs
        The Last Closed Flux Surfaces. Only used for bounds checking.

    Example
    -------
    ``` py
    >>> lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> qfactor = dex.UnityQfactor(lcfs)

    ```
    """

    _r: _PyQfactor

    def __init__(self, lcfs: LastClosedFluxSurface) -> None:
        self._r = _PyQfactor.build_unity(lcfs._r)
        super(EquilibriumObject, self).__init__()
        super(FluxCommute, self).__init__()
        super(Qfactor, self).__init__()

    def __repr__(self) -> str:
        return self._r.__repr__()

    def __str__(self) -> str:
        return self._r.__repr__()
