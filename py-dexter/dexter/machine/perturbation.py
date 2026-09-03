"""Defines the Perturbation wrapper over `_PyPerturbation`"""

import numpy as np
from collections.abc import Collection
from typing import TypeAlias

from .modes import FluteMode, NcFluteMode
from dexter._utils import _ReprStrImpl
from dexter._core import _PyPerturbation
from dexter.types import ArrayLike, Array

ModeObject: TypeAlias = FluteMode | NcFluteMode


class Perturbation(_ReprStrImpl):
    """A container type for [`ModeObjects`](modes), representing the total perturbation of the system.

    A `Perturbation` can consist of an arbitrary amount of modes, not necessarily of the same type.

    Parameters
    ----------
    modes
        Collection containing the comprising modes. If `None`, the perturbation is 0.

    Example
    -------
    ``` py title="Perturbation consisting of analytical flute modes"
    >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> per = dex.Perturbation(
    ...     [
    ...         dex.FluteMode(1e-4, LCFS, 3, 2, 0),
    ...         dex.FluteMode(2e-4, LCFS, 4, 3, 0),
    ...         dex.FluteMode(3e-4, LCFS, 5, 4, 0),
    ...     ]
    ... )

    ```

    Example
    -------
    ``` py title="Perturbation consisting of flute modes from a netCDF file"
    >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> per = dex.Perturbation(
    ...     [
    ...         dex.NcFluteMode("./netcdf.nc", "Cubic", 2, 1),
    ...         dex.NcFluteMode("./netcdf.nc", "Cubic", 3, 2),
    ...     ]
    ... )

    ```

    Example
    -------
    ``` py title="Perturbation consisting of mixed types of flute modes"
    >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> per = dex.Perturbation(
    ...     [
    ...         dex.FluteMode(1e-4, LCFS, 3, 2, 0),
    ...         dex.FluteMode(2e-4, LCFS, 5, 2, 0),
    ...         dex.NcFluteMode("./netcdf.nc", "Cubic", 3, 2),
    ...     ]
    ... )

    ```

    """

    _r: _PyPerturbation

    def __init__(self, modes: Collection[ModeObject] | None = None) -> None:
        if modes is not None:
            _modes = [mode._r for mode in modes]
        else:
            _modes = []
        self._r = _PyPerturbation(_modes)
        self._p_of_psi = np.vectorize(self._r.p_of_psi)
        self._p_of_psip = np.vectorize(self._r.p_of_psip)
        self._dp_dpsi = np.vectorize(self._r.dp_dpsi)
        self._dp_dpsip = np.vectorize(self._r.dp_dpsip)
        self._dp_of_psi_dtheta = np.vectorize(self._r.dp_of_psi_dtheta)
        self._dp_of_psip_dtheta = np.vectorize(self._r.dp_of_psip_dtheta)
        self._dp_of_psi_dzeta = np.vectorize(self._r.dp_of_psi_dzeta)
        self._dp_of_psip_dzeta = np.vectorize(self._r.dp_of_psip_dzeta)
        self._dp_of_psi_dt = np.vectorize(self._r.dp_of_psi_dt)
        self._dp_of_psip_dt = np.vectorize(self._r.dp_of_psip_dt)

    def p_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's value $p(\psi, \theta, \zeta, t)$, in Normalized Units."""
        return self._p_of_psi(psi, theta, zeta, t)[()]

    def p_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's value $p(\psi_p, \theta, \zeta, t)$, in Normalized Units."""
        return self._p_of_psip(psip, theta, zeta, t)[()]

    def dp_dpsi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi, \theta, \zeta, t)/d\psi$, in Normalized Units."""
        return self._dp_dpsi(psi, theta, zeta, t)[()]

    def dp_dpsip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi_p, \theta, \zeta, t)/d\psi_p$, in Normalized Units."""
        return self._dp_dpsip(psip, theta, zeta, t)[()]

    def dp_of_psi_dtheta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi, \theta, \zeta, t)/d\theta$, in Normalized Units."""
        return self._dp_of_psi_dtheta(psi, theta, zeta, t)[()]

    def dp_of_psip_dtheta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi_p, \theta, \zeta, t)/d\theta$, in Normalized Units."""
        return self._dp_of_psip_dtheta(psip, theta, zeta, t)[()]

    def dp_of_psi_dzeta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi, \theta, \zeta, t)/d\zeta$, in Normalized Units."""
        return self._dp_of_psi_dzeta(psi, theta, zeta, t)[()]

    def dp_of_psip_dzeta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi_p, \theta, \zeta, t)/d\zeta$, in Normalized Units."""
        return self._dp_of_psip_dzeta(psip, theta, zeta, t)[()]

    def dp_of_psi_dt(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi, \theta, \zeta, t)/dt$, in Normalized Units."""
        return self._dp_of_psi_dt(psi, theta, zeta, t)[()]

    def dp_of_psip_dt(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the perturbation's derivative $dp(\psi_p, \theta, \zeta, t)/dt$, in Normalized Units."""
        return self._dp_of_psip_dt(psip, theta, zeta, t)[()]

    def __len__(self) -> int:
        """Returns the total number of modes."""
        return self._r.__len__()
