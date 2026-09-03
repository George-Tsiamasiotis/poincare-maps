from typing import TypeAlias

from dexter.types import (
    Array,
    Array1,
    Array2,
    ArrayLike,
    ArrayShape,
    MachineType,
    NetCDFVersion,
    FluxCoordinate,
    FluxCoordinateState,
    Interpolation1dType,
    Interpolation2dType,
    PhaseMethod,
    CoordinateSet,
)


from dexter.machine.utils import LastClosedFluxSurface

from dexter.machine.geometries import LarGeometry, NcGeometry
from dexter.machine.qfactors import UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.machine.currents import LarCurrent, NcCurrent
from dexter.machine.bfields import LarBfield, NcBfield
from dexter.machine.modes import FluteMode, NcFluteMode
from dexter.machine.perturbation import Perturbation

from dexter.simulate.initial import InitialFlux, InitialConditions

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
    "MachineType",
    "NetCDFVersion",
    "FluxCoordinate",
    "FluxCoordinateState",
    "Interpolation1dType",
    "Interpolation2dType",
    "PhaseMethod",
    "CoordinateSet",
    # Machine
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
    # Simulate
    "InitialFlux",
    "InitialConditions",
]
