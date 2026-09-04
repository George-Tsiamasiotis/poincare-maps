"""Defines the different geometry objects as wrappers over `_PyGeometry`."""

import numpy as np
from semver import Version
from typing import TypeAlias

from dexter._core import _PyGeometry

from .utils import LastClosedFluxSurface
from .base import MachineObject, Geometry
from dexter.types import (
    Array1,
    Array2,
    ArrayShape,
    FluxCoordinateState,
    Interpolation1dType,
    Interpolation2dType,
    NetCDFVersion,
)


class LarGeometry(MachineObject, Geometry):
    r"""Analytical Large Aspect Ratio Geometry of a circular device.

    Parameters
    ----------
    baxis
        The magnetic field strength on the axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
    rlast
        The $r$ coordinate's value at the last closed flux surface in $[m]$.

    Notes
    -----
    + No $\psi/\psi_p$ bounds checks are performed in evaluations.

    + Evaluations methods that calculate or accept $\psi_p$ as a parameter always raise an Exception, since
    $\psi_p$ is defined through the q-factor.

    + The Jacobian is not available, since it is defined through $q$, $g$, $I$ and $B$.

    + In LAR equilibria, it holds that $R_0 \equiv R_{geo}$.

    + The definitions are not very strict at the moment.

    Example
    -------
    ``` py
    >>> geometry = dex.LarGeometry(
    ...     baxis=2, # Tesla
    ...     raxis=1.75, # meters
    ...     rlast=0.5, # meters
    ... )

    ```
    """

    _r: _PyGeometry

    def __init__(self, baxis: float, raxis: float, rlast: float) -> None:
        self._r = _PyGeometry.build_lar(baxis, raxis, rlast)
        super(MachineObject, self).__init__()
        super(Geometry, self).__init__()


class NcGeometry(MachineObject, Geometry):
    r"""Geometry of a realistic configuration.

    Stores fluxes, angles and lab variables’ data, and provides interpolation methods between them.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp1d_type
        The 1D interpolation type for the 1D quantities.
    interp2d_type
        The 2D interpolation type for the 2D quantities.

    Example
    -------
    ``` py
    >>> geometry = dex.NcGeometry(path, "Cubic", "Bicubic")

    ```
    """

    _r: _PyGeometry

    def __init__(
        self,
        path: str,
        interp1d_type: Interpolation1dType,
        interp2d_type: Interpolation2dType,
    ) -> None:
        self._r = _PyGeometry.build_nc(path, interp1d_type, interp2d_type)
        super(MachineObject, self).__init__()
        super(Geometry, self).__init__()

    @property
    def path(self) -> str:
        """The path to the NetCDF file."""
        return self._r.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The path to the NetCDF file."""
        return Version.parse(self._r.netcdf_version)

    @property
    def interp1d_type(self) -> Interpolation1dType:
        """The 1D interpolation type."""
        return self._r.interp1d_type

    @property
    def interp2d_type(self) -> Interpolation2dType:
        """The 2D interpolation type."""
        return self._r.interp2d_type

    @property
    def shape(self) -> ArrayShape:
        r"""Returns the $(\psi/\psi_p,\theta)$ shape of the 2D arrays."""
        return self._r.shape

    @property
    def psi_array(self) -> Array1:
        """The toroidal flux's values."""
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        """The poloidal flux's values."""
        return self._r.get_array("psip_array")

    @property
    def theta_array(self) -> Array1:
        r"""The poloidal angle $\theta$ values."""
        return self._r.get_array("theta_array")

    @property
    def r_array(self) -> Array1:
        r"""The radial coordinate $r$ values, in $[m]$."""
        return self._r.get_array("r_array")

    @property
    def rlab_array(self) -> Array2:
        r"""The $R$ values, in $[m]$."""
        return self._r.get_array2d("rlab_array")

    @property
    def zlab_array(self) -> Array2:
        r"""The $Z$ values, in $[m]$."""
        return self._r.get_array2d("zlab_array")

    @property
    def jacobian_array(self) -> Array2:
        r"""The Jacobian $J$ values, in $[m]$."""
        return self._r.get_array2d("jacobian_array")


GeometryObject: TypeAlias = LarGeometry | NcGeometry
r"""Available [`Geometry`][dexter.machine.base.Geometry] objects.

+ [`LarGeometry`][dexter.LarGeometry]: Analytical Large Aspect Ratio Geometry of a circular device.
+ [`NcGeometry`][dexter.NcGeometry]: Geometry of a realistic configuration.
"""
