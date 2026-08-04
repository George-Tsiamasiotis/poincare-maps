"""Defines the different plasma current objects as wrappers over `_PyCurrent`."""

import numpy as np
from semver import Version

from dexter._core import _PyCurrent

from .utils import LastClosedFluxSurface
from .base import EquilibriumObject, Current
from dexter.types import Array1, FluxCoordinateState, Interpolation1dType, NetCDFVersion


class LarCurrent(EquilibriumObject, Current):
    r"""Analytical Large Aspect Ratio Current with $g=1$ and $I=0$.

    !!! note

        No $\psi/\psi_p$ bounds checks are performed in evaluations.

    Example
    -------
    ``` py
    >>> current = dex.LarCurrent()

    ```
    """

    _r: _PyCurrent

    def __init__(self) -> None:
        self._r = _PyCurrent.build_lar()
        super(EquilibriumObject, self).__init__()
        super(Current, self).__init__()


class NcCurrent(EquilibriumObject, Current):
    r"""Numerical plasma current profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D interpolation type.

    Example
    -------
    ``` py
    >>> current = dex.NcCurrent(path, "Cubic")

    ```
    """

    _r: _PyCurrent

    def __init__(self, path: str, interp_type: Interpolation1dType) -> None:
        self._r = _PyCurrent.build_nc(path, interp_type)
        super(EquilibriumObject, self).__init__()
        super(Current, self).__init__()

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
    def psi_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{LCFS}$."""
        return self._r.psi_last

    @property
    def psip_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{p,LCFS}$."""
        return self._r.psip_last

    @property
    def psi_array(self) -> Array1:
        """The toroidal flux's values."""
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        """The poloidal flux's values."""
        return self._r.get_array("psip_array")

    @property
    def g_array(self) -> Array1:
        """The poloidal plasma current $g$ values."""
        return self._r.get_array("g_array")

    @property
    def i_array(self) -> Array1:
        """The toroidal plasma current $I$ values."""
        return self._r.get_array("i_array")
