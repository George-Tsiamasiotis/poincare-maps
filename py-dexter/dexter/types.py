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

# =================== Equilibrium

ObjectType: TypeAlias = Literal["Analytical", "Numerical"]
""" Describes the type of equilibrium the object represents.

# Numerical
Equilibrium reconstructed from numerical data. Evaluations are calculated by interpolation.

# Analytical
Equilibrium described by analytical formulas. Evaluations are calculated by simply evaluating
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
