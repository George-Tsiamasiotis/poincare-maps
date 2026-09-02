from dexter._core import _PyLastClosedFluxSurface

from dexter.types import FluxCoordinate
from dexter._utils import _ReprStrImpl


class LastClosedFluxSurface(_ReprStrImpl):
    """Helper type to define the Last Closed Flux Surface (LCFS) with respect to one of the two fluxes.

    This type is instantiated through the [`Toroidal`][dexter.LastClosedFluxSurface.Toroidal] and
    [`Poloidal`][dexter.LastClosedFluxSurface.Poloidal] class methods.
    """

    _r: _PyLastClosedFluxSurface

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def Toroidal(cls, value: float) -> LastClosedFluxSurface:
        r"""Defines the Last Closed Flux Surface with respect to $\psi$.

        Parameters
        ----------
        value
            The value of the toroidal magnetic flux at the last closed flux surface.

        Example
        -------
        ```python title="LastClosedFluxSurface creation"
        >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)

        ```
        """
        obj = LastClosedFluxSurface.__new__(LastClosedFluxSurface)
        obj._r = _PyLastClosedFluxSurface.toroidal(value)
        return obj

    @classmethod
    def Poloidal(cls, value: float) -> LastClosedFluxSurface:
        r"""Defines the Last Closed Flux Surface with respect to $\psi$.

        Parameters
        ----------
        value
            The value of the toroidal magnetic flux at the last closed flux surface.

        Example
        -------
        ```python title="LastClosedFluxSurface creation"
        >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.05)

        ```
        """
        obj = LastClosedFluxSurface.__new__(LastClosedFluxSurface)
        obj._r = _PyLastClosedFluxSurface.poloidal(value)
        return obj

    @property
    def value(self) -> float:
        """The value of the magnetic flux."""
        return self._r.value

    @property
    def kind(self) -> FluxCoordinate:
        """The kind of the magnetic flux."""
        return self._r.kind

    @classmethod
    def _wrap(cls, _r: _PyLastClosedFluxSurface) -> LastClosedFluxSurface:
        new = _PyLastClosedFluxSurface.__new__(cls)
        new._r = _r
        return new
