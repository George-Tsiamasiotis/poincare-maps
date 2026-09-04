r"""Type Aliases used throughout the package."""

import numpy as np
from typing import TypeAlias, Literal
from collections.abc import Sequence

from semver import Version

# =================== Common


ArrayShape: TypeAlias = tuple[int, ...]
"""Shape of a numpy array."""

Array: TypeAlias = np.ndarray[ArrayShape, np.dtype[np.float64]]
"""ND numpy array."""

Array1: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
"""1D numpy array."""

Array2: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]
"""2D numpy array."""

ArrayLike: TypeAlias = float | Array | Sequence
"""Objects that can be converted to arrays, i.e. float, np.ndarray, sequences, ..."""

ParticleSpecies: TypeAlias = Literal[
    "Electron",
    "Proton",
    "Deuterium",
    "Tritium",
    "Alpha",
    "He3",
]
""" The particle species under study. This only affects unit conversions to SI."""

Unit: TypeAlias = (
    Literal[
        "NormMeter",
        "NormTesla",
        "NormSecond",
        "NormHertz",
        "NormJoule",
    ]
    | str
)
"""Strings parsed by `pint` as units, with the normalized units added."""

# =================== Machine

MachineType: TypeAlias = Literal["Analytical", "Numerical"]
""" Describes the type of machine the object represents.

# Numerical
Configuration reconstructed from numerical data. Evaluations are calculated by interpolation.

# Analytical
Machine described by analytical formulas. Evaluations are calculated by simply evaluating
the formulas.
"""

NetCDFVersion: TypeAlias = Version
"""The netCDF convention version (SemVer)."""

FluxCoordinate: TypeAlias = Literal["Toroidal", "Poloidal"]
r"""Magnetic flux coordinates $\psi$ and $\psi_p$."""

FluxCoordinateState: TypeAlias = Literal["Good", "Bad", "None"]
"""Flux coordinate state.

# Good
Values exist and are increasing. Can be used as an x coordinate for evaluations.

# Bad
Values exist but are not monotonic. Can only be used as y-data.

#None
Variable does not exist or is empty.
"""

Interpolation1dType: TypeAlias = Literal[
    "Linear", "Cubic", "Cubic Periodic", "Akima", "Akima Periodic", "Steffen"
]
"""Available 1D Interpolation types

# Linear
Simple linear interpolation

# Cubic
Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on each
interval, with matching first and second derivatives at the supplied data-points. The second
derivative is chosen to be zero at the first point and last point.

# CubicPeriodic
Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on each
interval, with matching first and second derivatives at the supplied data-points. The derivatives
at the first and last points are also matched. Note that the last point in the data must have the
same y-value as the first point, otherwise the resulting periodic interpolation will have a
discontinuity at the boundary.

# Akima
Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded corner
algorithm of Wodicka.

# AkimaPeriodic
Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded corner
algorithm of Wodicka.

# Steffen
Steffen interpolation. Steffen’s method guarantees the monotonicity of the interpolating function
between the given data points. Therefore, minima and maxima can only occur exactly at the data
points, and there can never be spurious oscillations between data points. The interpolated function
is piecewise cubic in each interval. The resulting curve and its first derivative are guaranteed
to be continuous, but the second derivative may be discontinuous.
"""

Interpolation2dType: TypeAlias = Literal["Bilinear", "Bicubic"]
"""Available 2D Interpolation Types.

# Bilinear
Simple bilinear interpolation.

# Bicubic
Bicubic Interpolation.
"""

PhaseMethod: TypeAlias = (
    Literal["Zero", "Average", "Resonance", "Interpolation"]
    | tuple[Literal["Custom"], float]
)
r""" Defines the calculation method of the phase $\phi$ in a Numerical Harmonic.

# Zero
Corresponds to $\phi = 0$.

# Average
Corresponds to $\phi = const =$ the average of all the values of the `phase_array`.

# Resonance
Corresponds to $\phi = const =$ the value of $\phi$ at the resonance $m/n$. In the case that the
resonance falls outside the last closed flux surface, or does not correspond to a valid q-factor
value, it defaults to `Zero`.

# Interpolation
Interpolation over the `phase_array`.

# Custom(float)
Use a custom value for $\phi = const$.
"""

# =================== Simulate

CoordinateSet: TypeAlias = Literal[
    "BoozerToroidal",
    "BoozerPoloidal",
    "MixedToroidal",
    "MixedPoloidal",
]
r""" The kind of InitialConditions set.

    - `BoozerToroidal`: Initial conditions set in the $(t, \psi, \theta, \zeta, \rho, \mu)$ space.
    - `BoozerPoloidal`: Initial conditions set in the $(t, \psi_p, \theta, \zeta, \rho, \mu)$ space.
    - `MixedToroidal`: Initial conditions set in the $(t, P_\zeta, \psi, \theta, \zeta, \mu)$ space.
    - `MixedPoloidal`: Initial conditions set in the $(t, P_\zeta, \psi_p, \theta, \zeta, \mu)$ space.
"""

SteppingMethod = (
    Literal["EnergyAdaptiveStep", "ErrorAdaptiveStep"]
    | tuple[Literal["FixedStep"], float]
)
"""The stepping method of the solver.

- `EnergyAdaptiveStep`: Forces the step size to be small enough so that the Energy difference
    from step to step is under a certain threshold. The tolerances can be adjusted with the
    energy_rel_tol and energy_abs_tol fields.
- `ErrorAdaptiveStep`: Classic RK error estimation : Adjust the step size to minimize the
    local truncation error.
- `FixedStep(float)`: Fixed step size.
"""


Intersection: TypeAlias = Literal["ConstZeta", "ConstTheta"]
r""" Defines the surface of the Poincare section.

- `ConstTheta`: Defines a surface of $\theta = const$.
- `ConstZeta`: Defines a surface of $\zeta = const$.
"""

IntegrationStatus: TypeAlias = Literal[
    "Initialized",
    "PartlyInitialized",
    "InvalidInitialConditions",
    "OutOfBoundsInitialization",
    "Integrated",
    "Intersected",
    "ClosedPeriods(..)",
    "Escaped",
    "ModStateEscaped",
    "IntersectedTimedOut",
    "InvalidIntersections",
    "TimedOut(...)",
    "Failed(...)",
]
r"""The integration status of a Particle.

- `Initialized`: Initialized by InitialConditions, not integrated.
- `PartlyInitialized`: InitialConditions have not been fully calculated yet.
- `InvalidInitialConditions`: Invalid InitialConditions. May occur when using Mixed variables
    with objects that cannot define them, for example Mixed Toroidal coordinates when $g(\psi)$ is
    not defined.
- `OutOfBoundsInitialization`: InitialConditions where out of bounds.
- `Integrated`: Reached the end of the integration successfully.
- `Intersected`: Intersections calculation successful.
- `ClosedPeriods(..)`: Integrated for a certain amount of θ-ψ periods.
- `Escaped`: Escaped the last closed flux surface (LCFS).
- `ModStateEscaped`: Escaped when performing a step on the modified system. This indicates that
    something is wrong in Hénon’s trick implementation.
- `IntersectedTimedOut`: Calculated some intersections correctly but also timed out.
- `InvalidIntersections`: Calculated invalid intersections.
- `TimedOut(...)`: Timed out after a maximum number of steps.
- `Failed(...)`: Simulation failed for unknown reasons.
"""

EnergyPzetaPosition = Literal[
    "Alpha",
    "Beta",
    "Gamma",
    "Delta",
    "Epsilon",
    "Zeta",
    "Eta",
    "Theta",
    "Iota",
    "Kappa",
    "Lambda",
    "Mu",
    "Unclassified",
]
r"""The position of an $(E, P_\zeta)$ point on the $(E, P_\zeta)$ plane, relative to the orbit
classification curves.

See the diagram for explanation.
"""

OrbitType: TypeAlias = Literal[
    "Undefined",
    "TrappedLost",
    "TrappedConfined",
    "CoPassingLost",
    "CoPassingConfined",
    "CuPassingLost",
    "CuPassingConfined",
    "Potato",
    "Stagnated",
    "Unclassified",
    "Failed(..)",
]
r"""A particle's orbit type, calculated through the [`dexter.Particle.close()`] routine.

- `Undefined`: Particle has not been classified.
- `TrappedLost`: A Trapped-Lost particle. A particle is called trapped if there exists a
    mirror point where $\rho=0$.
- `TrappedConfined`: A Trapped-Confined particle. A particle is called trapped if there
    exists a mirror point where $\rho=0$.
- `CoPassingLost`: A CoPassing-Lost particle. A particle is called co-passing if it is not
    trapped and it holds that $\dot\theta>0$.
- `CoPassingConfined`: A CoPassing-Confined particle. A particle is called co-passing if it
    is not trapped and it holds that $\dot\theta>0$.
- `CuPassingLost`: A CounterPassing-Lost particle. A particle is called counter-passing if it is
    not trapped and it holds that $\dot\theta<0$.
- `CuPassingConfined`: A CounterPassing-Confined particle. A particle is called counter-passing
    if it is not trapped and it holds that $\dot\theta<0$.
- `Potato`: A Potato particle. A particle’s orbit is called a potato orbit if it is trapped
    but still circles the magnetic axis due to its drift. In the $(E, P_\zeta)$ plane, those
    lie inside the intersection of the trapped-passing boundary and the magnetic axis parabola.
- `Stagnated`: A Potato particle. A particle is called stagnated if it always has positive
    parallel velocity but does not circle the magnetic axis. In the $(E, P_\zeta)$ plane,
    those lie to the right of the trapped-passing boundary and above the magnetic axis parabola.
- `Unclassified`: Not falling under any of the other categories.
- `Failed(..)`: Error classifying the orbit.
"""
