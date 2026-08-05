"""Defines the different Mode objects as wrappers over `_PyMode`."""

import numpy as np
from semver import Version

from dexter._core import _PyMode

from .utils import LastClosedFluxSurface
from .base import EquilibriumObject, Mode
from dexter.types import (
    Array1,
    FluxCoordinateState,
    Interpolation1dType,
    NetCDFVersion,
    PhaseMethod,
)


class FluteMode(EquilibriumObject, Mode):
    r"""A simple analytical flute mode.

    Parameters
    ----------
    epsilon
        The modes's "amplitude" $\epsilon$. Corresponds the value of the amplitude at the last
        closed flux surface.
    lcfs
        The Last Closed Flux Surface, with respect to which the mode is defined.
    m
        The poloidal mode number.
    n
        The toroidal mode number.
    phase
        The mode's constant phase $\phi$.

    Example
    -------
    ``` py
    >>> lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> mode = dex.FluteMode(1e-4, lcfs, 3, 4, 0)

    ```
    """

    _r: _PyMode

    def __init__(
        self, epsilon: float, lcfs: LastClosedFluxSurface, m: int, n: int, phase: float
    ) -> None:
        self._r = _PyMode.build_flute(epsilon, lcfs._r, m, n, phase)
        super(EquilibriumObject, self).__init__()
        super(Mode, self).__init__()

    @property
    def lcfs(self) -> LastClosedFluxSurface:
        r"""The mode's Last Closed Flux Surface."""
        return LastClosedFluxSurface._wrap(self._r.lcfs)

    @property
    def epsilon(self) -> float:
        r"""The mode's constant 'amplitude' $\epsilon$."""
        return self._r.epsilon

    @property
    def phase(self) -> float:
        r"""The mode's constant phase $\phi$."""
        return self._r.phase


class NcFluteMode(EquilibriumObject, Mode):
    r"""Single perturbation flute mode from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D interpolation type.
    m
        The poloidal mode number.
    n
        The toroidal mode number.
    phase_method
        The phase $\phi(\psi/\psi_p)$ calculation method.
    analytical_threshold_index
        The modes’s analytical threshold point. Defines the index of the magnetic flux values under
        which to switch to the analytica formula.

        !!! note "Numerical flute mode analytical patching"

            By definition, flute modes must behave like $\approx\sqrt\psi$ close to the axis, and
            therefore their derivative with respect to the flux must go to infinity. This is a
            behavior that splines cannot replicate, resulting to unnatural orbits close to the
            magnetic axis. To solve this, the mode switches to an analytical formula for the values
            of $\psi/\psi_p$ under a certain threshold. The threshold is defined by the flux value
            at the position index of the data array.

            #Formula
            The patch has the form $\beta\sqrt\psi + \gamma$, where $\beta$ and $\gamma$ are
            adjusted in order to ensure continuity of both $\alpha(\psi)$ and its first derivative.
            $\beta$ is calculated first by $\beta = 2\alpha' \sqrt\psi$ to ensure the correct
            value of the derivative $d\alpha/d\psi$ at the patch’s edge. Finally,
            $\gamma = \alpha - \beta\sqrt\psi$ ensures the continuity of $\alpha$ itself.

            Note that sometimes $\gamma$ may become slightly negative, resulting to $\alpha$
            becoming slightly negative extremely close to the axis. However this error should be
            negligible compared to the possible non-continuity of $\alpha$’s higher derivatives
            and/or its deviation from the actual data.

    Example
    -------
    ``` py
    >>> mode = dex.NcFluteMode(path, "Cubic", 3, 2, phase_method="Zero")
    >>> mode = dex.NcFluteMode(
    ...     path=path,
    ...     interp_type="Cubic",
    ...     m=3,
    ...     n=2,
    ...     phase_method=("Custom", 1.57),
    ...     analytical_threshold_index=0,
    ... )

    ```
    """

    _r: _PyMode

    def __init__(
        self,
        path: str,
        interp_type: Interpolation1dType,
        m: int,
        n: int,
        phase_method: PhaseMethod = "Interpolation",
        analytical_threshold_index: int = 3,
    ) -> None:
        self._r = _PyMode.build_nc(
            path, interp_type, m, n, phase_method, analytical_threshold_index
        )
        super(EquilibriumObject, self).__init__()
        super(Mode, self).__init__()

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
    def phase_method(self) -> PhaseMethod:
        r"""The phase $\phi$ calculation method."""
        return self._r.phase_method

    @property
    def analytical_threshold_index(self) -> int:
        r"""The analytical threshold index."""
        return self._r.analytical_threshold_index

    @property
    def phase_average(self) -> float | None:
        r"""The average value of the phase array, if `PhaseMethod == 'Average'`."""
        return self._r.phase_average if self.phase_method == "Average" else None

    @property
    def psi_array(self) -> Array1:
        """The toroidal flux's values."""
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        """The poloidal flux's values."""
        return self._r.get_array("psip_array")

    @property
    def alpha_array(self) -> Array1:
        r"""The $\alpha$ values."""
        return self._r.get_array("alpha_array")

    @property
    def phase_array(self) -> Array1:
        r"""The $\phase$ values."""
        return self._r.get_array("phase_array")
