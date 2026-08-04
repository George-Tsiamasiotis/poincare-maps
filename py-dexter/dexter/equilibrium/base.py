"""Equilibrium objects' base classes.

Equilibrium objects define evaluations over equilibrium quantities, provide information about the
[`state`][dexter.types.FluxCoordinateState] of each magnetic flux coordinate, as well as useful
scalar quantities and data arrays.

Each parent class corresponds to an evaluation Trait on the Rust API.
"""

import numpy as np

from dexter._core import _PyQfactor, _PyCurrent
from dexter._utils import _ReprStrImpl
from dexter.types import ArrayLike, Array, FluxCoordinateState, ObjectType


class EquilibriumObject(_ReprStrImpl):
    """Common attributes in all equilibrium objects."""

    _r: _PyQfactor | _PyCurrent

    @property
    def object_type(self) -> ObjectType:
        """The object’s equilibrium type."""
        return self._r.object_type

    @property
    def psi_state(self) -> FluxCoordinateState:
        r"""The state of the toroidal flux coordinate $\psi$."""
        return self._r.psi_state

    @property
    def psip_state(self) -> FluxCoordinateState:
        r"""The state of the toroidal flux coordinate $\psi_p$."""
        return self._r.psip_state


class FluxCommute(_ReprStrImpl):
    """Methods for converting from one magnetic flux to the other."""

    _r: _PyQfactor

    def __init__(self) -> None:
        self._psi_of_psip = np.vectorize(self._r.psi_of_psip)
        self._psip_of_psi = np.vectorize(self._r.psip_of_psi)

    def psip_of_psi(self, psi: ArrayLike) -> Array:
        r"""Calculates $\psi_p(\psi)$, in Normalized Units."""
        return self._psip_of_psi(psi)[()]

    def psi_of_psip(self, psip: ArrayLike) -> Array:
        r"""Calculates $\psi(\psi_p)$, in Normalized Units."""
        return self._psi_of_psip(psip)[()]


class Qfactor(_ReprStrImpl):
    """q-factor related quantities and evaluation methods."""

    _r: _PyQfactor

    def __init__(self) -> None:
        self._q_of_psi = np.vectorize(self._r.q_of_psi)
        self._q_of_psip = np.vectorize(self._r.q_of_psip)
        self._dpsip_dpsi = np.vectorize(self._r.dpsip_dpsi)
        self._dpsi_dpsip = np.vectorize(self._r.dpsi_dpsip)
        self._psi_of_q = np.vectorize(self._r.psi_of_q)
        self._psip_of_q = np.vectorize(self._r.psip_of_q)
        self._iota_of_psi = np.vectorize(self._r.iota_of_psi)
        self._iota_of_psip = np.vectorize(self._r.iota_of_psip)

    @property
    def psi_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{LCFS}$."""
        return self._r.psi_last

    @property
    def psip_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{p,LCFS}$."""
        return self._r.psip_last

    @property
    def qlast(self) -> float:
        r"""The q-factor's value at the last closed flux surface, $q_{LCFS}$."""
        return self._r.qlast

    @property
    def qaxis(self) -> float:
        r"""The q-factor's value at the magnetic axis, $q_{axis}$."""
        return self._r.qaxis

    def q_of_psi(self, psi: ArrayLike) -> Array:
        r"""Calculates $q(\psi)$, in Normalized Units."""
        return self._q_of_psi(psi)[()]

    def q_of_psip(self, psip: ArrayLike) -> Array:
        r"""Calculates $q(\psi_p)$, in Normalized Units."""
        return self._q_of_psi(psip)[()]

    def dpsip_dpsi(self, psi: ArrayLike) -> Array:
        r"""Calculates $d\psi_p/d\psi$, in Normalized Units."""
        return self._dpsip_dpsi(psi)[()]

    def dpsi_dpsip(self, psip: ArrayLike) -> Array:
        r"""Calculates $d\psi/d\psi_p$, in Normalized Units."""
        return self._dpsi_dpsip(psip)[()]

    def psi_of_q(self, q: ArrayLike) -> Array:
        r"""Calculates $\psi(q)$, in Normalized Units."""
        return self._psi_of_q(q)[()]

    def psip_of_q(self, q: ArrayLike) -> Array:
        r"""Calculates $\psi_p(q)$, in Normalized Units."""
        return self._psip_of_q(q)[()]

    def iota_of_psi(self, psi: ArrayLike) -> Array:
        r"""Calculates $\iota(\psi)$, in Normalized Units."""
        return self._iota_of_psi(psi)[()]

    def iota_of_psip(self, psip: ArrayLike) -> Array:
        r"""Calculates $\iota(\psi_p)$, in Normalized Units."""
        return self._iota_of_psip(psip)[()]


class Current(_ReprStrImpl):
    """Plasma current related evaluation methods."""

    _r: _PyCurrent

    def __init__(self) -> None:
        self._g_of_psi = np.vectorize(self._r.g_of_psi)
        self._g_of_psip = np.vectorize(self._r.g_of_psip)
        self._i_of_psi = np.vectorize(self._r.i_of_psi)
        self._i_of_psip = np.vectorize(self._r.i_of_psip)
        self._dg_dpsi = np.vectorize(self._r.dg_dpsi)
        self._dg_dpsip = np.vectorize(self._r.dg_dpsip)
        self._di_dpsi = np.vectorize(self._r.di_dpsi)
        self._di_dpsip = np.vectorize(self._r.di_dpsip)

    def g_of_psi(self, psi: ArrayLike) -> Array:
        r"""Calculates $g(\psi)$, in Normalized Units."""
        return self._g_of_psi(psi)[()]

    def g_of_psip(self, psip: ArrayLike) -> Array:
        r"""Calculates $g(\psi_p)$, in Normalized Units."""
        return self._g_of_psip(psip)[()]

    def i_of_psi(self, psi: ArrayLike) -> Array:
        r"""Calculates $I(\psi)$, in Normalized Units."""
        return self._i_of_psi(psi)[()]

    def i_of_psip(self, psip: ArrayLike) -> Array:
        r"""Calculates $I(\psi_p)$, in Normalized Units."""
        return self._i_of_psip(psip)[()]

    def dg_dpsi(self, psi: ArrayLike) -> Array:
        r"""Calculates $dg/d\psi$, in Normalized Units."""
        return self._dg_dpsi(psi)[()]

    def dg_dpsip(self, psip: ArrayLike) -> Array:
        r"""Calculates $dg/d\psi_p$, in Normalized Units."""
        return self._dg_dpsip(psip)[()]

    def di_dpsi(self, psi: ArrayLike) -> Array:
        r"""Calculates $dI/d\psi$, in Normalized Units."""
        return self._di_dpsi(psi)[()]

    def di_dpsip(self, psip: ArrayLike) -> Array:
        r"""Calculates $dI/d\psi_p$, in Normalized Units."""
        return self._di_dpsip(psip)[()]
