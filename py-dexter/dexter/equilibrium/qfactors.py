"""Defines the different q-factor objects as wrappers over `_PyQfactor`."""

import numpy as np
from semver import Version

from dexter._core import _PyQfactor

from .utils import LastClosedFluxSurface
from .base import EquilibriumObject, FluxCommute, Qfactor
from dexter.types import Array1, FluxCoordinateState, Interpolation1dType, NetCDFVersion


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


class ParabolicQfactor(EquilibriumObject, FluxCommute, Qfactor):
    r"""Analytical parabolic q-factor profile.

    Parameters
    ----------
    qaxis
        The q-factor's value at the magnetic axis, $q_{axis}$.
    qwall
        The q-factor's value at the last closed flux surface, $q_{LCFS}$.

    Example
    -------
    ``` py
    >>> lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> qfactor = dex.ParabolicQfactor(1.1, 3.9, lcfs)

    ```
    """

    _r: _PyQfactor

    def __init__(self, qaxis: float, qwall: float, lcfs: LastClosedFluxSurface) -> None:
        self._r = _PyQfactor.build_parabolic(qaxis, qwall, lcfs._r)
        super(EquilibriumObject, self).__init__()
        super(FluxCommute, self).__init__()
        super(Qfactor, self).__init__()


class NcQfactor(EquilibriumObject, FluxCommute, Qfactor):
    r"""Numerical q-factor profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    If either psi_norm or psip_norm is missing from the netCDF file, it is calculated from the
    other by integrating $q(\psi_p)$ or $\iota(\psi)$ respectively. In the case that the calculated
    values are monotonic, the other flux can be used as a flux coordinate as well.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D interpolation type.

    Example
    -------
    ``` py
    >>> qfactor = dex.NcQfactor(path, "Cubic")

    ```
    """

    _r: _PyQfactor

    def __init__(self, path: str, interp_type: Interpolation1dType) -> None:
        self._r = _PyQfactor.build_nc(path, interp_type)
        super(EquilibriumObject, self).__init__()
        super(FluxCommute, self).__init__()
        super(Qfactor, self).__init__()

    @property
    def path(self) -> str:
        """The path to the NetCDF file."""
        return self._r.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The path to the NetCDF file."""
        return Version.parse(self._r.netcdf_version)

    @property
    def interp_type(self) -> Interpolation1dType:
        """The 1D interpolation type."""
        return self._r.interp_type

    @property
    def psi_array(self) -> Array1:
        """The toroidal flux's values."""
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        """The poloidal flux's values."""
        return self._r.get_array("psip_array")

    @property
    def q_array(self) -> Array1:
        """The q-factor's values."""
        return self._r.get_array("q_array")
