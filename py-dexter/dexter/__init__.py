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
    SteppingMethod,
    ParticleSpecies,
    IntegrationStatus,
    EnergyPzetaPosition,
    OrbitType,
)


from dexter.machine.utils import LastClosedFluxSurface

from dexter.machine.geometries import GeometryObject, LarGeometry, NcGeometry
from dexter.machine.qfactors import (
    QfactorObject,
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
)
from dexter.machine.currents import CurrentObject, LarCurrent, NcCurrent
from dexter.machine.bfields import BfieldObject, LarBfield, NcBfield
from dexter.machine.modes import ModeObject, FluteMode, NcFluteMode
from dexter.machine.perturbation import Perturbation
from dexter.machine.machine import Machine

from dexter.simulate.initial import InitialFlux, InitialConditions
from dexter.simulate.particle import Particle

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
    "SteppingMethod",
    "ParticleSpecies",
    "IntegrationStatus",
    "EnergyPzetaPosition",
    "OrbitType",
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
    "Machine",
    # Simulate
    "InitialFlux",
    "InitialConditions",
    "Particle",
]
