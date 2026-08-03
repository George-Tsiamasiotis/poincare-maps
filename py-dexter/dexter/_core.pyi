"""This file mirrors all the definitions made in the `py-dexter` Rust API."""

from dexter.types import FluxCoordinate

class _PyLastClosedFluxSurface:

    value: float
    kind: FluxCoordinate

    @classmethod
    def toroidal(cls, value: float) -> _PyLastClosedFluxSurface: ...
    @classmethod
    def poloidal(cls, value: float) -> _PyLastClosedFluxSurface: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
