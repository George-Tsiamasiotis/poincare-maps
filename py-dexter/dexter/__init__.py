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
    Interpolation2dType,
    PhaseMethod,
)


from dexter.equilibrium.utils import LastClosedFluxSurface

from dexter.equilibrium.geometries import LarGeometry, NcGeometry
from dexter.equilibrium.qfactors import UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.equilibrium.currents import LarCurrent, NcCurrent
from dexter.equilibrium.bfields import LarBfield, NcBfield
from dexter.equilibrium.modes import FluteMode, NcFluteMode
from dexter.equilibrium.perturbation import Perturbation

GeometryObject: TypeAlias = LarGeometry | NcGeometry
QfactorObject: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
CurrentObject: TypeAlias = LarCurrent | NcCurrent
BfieldObject: TypeAlias = LarBfield | NcBfield
ModeObject: TypeAlias = FluteMode | NcFluteMode

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
    "Interpolation2dType",
    "PhaseMethod",
    # Equilibrium
    "GeometryObject",
    "QfactorObject",
    "CurrentObject",
    "BfieldObject",
    "ModeObject",
    "LastClosedFluxSurface",
    "LarGeometry",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
    "LarBfield",
    "NcBfield",
    "FluteMode",
    "NcFluteMode",
    "Perturbation",
]
