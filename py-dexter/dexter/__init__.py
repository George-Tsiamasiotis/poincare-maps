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

QfactorObject: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor

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
    "LastClosedFluxSurface",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
]
