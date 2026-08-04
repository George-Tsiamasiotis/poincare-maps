from typing import TypeAlias

from dexter.types import (
    Array,
    Array1,
    Array2,
    ArrayLike,
    ArrayShape,
    ObjectType,
    NetCDFVersion,
    FluxCoordinate,
    FluxCoordinateState,
    Interpolation1dType,
)


from dexter.equilibrium.utils import LastClosedFluxSurface

from dexter.equilibrium.qfactors import UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.equilibrium.currents import LarCurrent, NcCurrent
from dexter.equilibrium.bfields import LarBfield, NcBfield

QfactorObject: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
CurrentObject: TypeAlias = LarCurrent | NcCurrent
BfieldObject: TypeAlias = LarBfield | NcBfield

__all__ = [
    # Type Aliases
    "Array",
    "Array1",
    "Array2",
    "ArrayLike",
    "ArrayShape",
    "ObjectType",
    "NetCDFVersion",
    "FluxCoordinate",
    "FluxCoordinateState",
    "Interpolation1dType",
    # Equilibrium
    "QfactorObject",
    "CurrentObject",
    "BfieldObject",
    "LastClosedFluxSurface",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
    "LarBfield",
    "NcBfield",
]
