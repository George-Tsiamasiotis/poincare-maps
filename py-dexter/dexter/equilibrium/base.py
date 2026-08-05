"""Equilibrium objects' base classes.

Equilibrium objects define evaluations over equilibrium quantities, provide information about the
[`state`][dexter.types.FluxCoordinateState] of each magnetic flux coordinate, as well as useful
scalar quantities and data arrays.

Each parent class corresponds to an evaluation Trait on the Rust API.

Classes
-------
EquilibriumObject
    Common attributes in all equilibrium objects.
FluxCommute
    Methods for converting from one magnetic flux to the other.
Qfactor
    q-factor related quantities and evaluation methods.
Current
    Plasma current related evaluation methods.
Bfield
    Magnetic field related evaluation methods.
Geometry
    Device geometry related evaluation methods.
Mode
    Single perturbation mode related evaluation methods.

"""

import numpy as np
from typing import Any

from dexter._core import _PyQfactor, _PyCurrent, _PyBfield, _PyGeometry, _PyMode
from dexter._utils import _ReprStrImpl
from dexter.types import ArrayLike, Array, Array1, FluxCoordinateState, ObjectType


class EquilibriumObject(_ReprStrImpl):
    """Common attributes in all equilibrium objects."""

    _r: Any

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


class Bfield(_ReprStrImpl):
    """Magnetic field related evaluation methods."""

    _r: _PyBfield

    def __init__(self) -> None:
        self._b_of_psi = np.vectorize(self._r.b_of_psi)
        self._b_of_psip = np.vectorize(self._r.b_of_psip)
        self._db_dpsi = np.vectorize(self._r.db_dpsi)
        self._db_dpsip = np.vectorize(self._r.db_dpsip)
        self._db_of_psi_dtheta = np.vectorize(self._r.db_of_psi_dtheta)
        self._db_of_psip_dtheta = np.vectorize(self._r.db_of_psip_dtheta)

    def b_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $B(\psi, \theta)$, in Normalized Units."""
        return self._b_of_psi(psi, theta)[()]

    def b_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $B(\psi_p, \theta)$, in Normalized Units."""
        return self._b_of_psip(psip, theta)[()]

    def db_dpsi(self, psi: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $dB(\psi, \theta)/d\psi$, in Normalized Units."""
        return self._db_dpsi(psi, theta)[()]

    def db_dpsip(self, psip: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $dB(\psi_p, \theta)/d\psi_p$, in Normalized Units."""
        return self._db_dpsip(psip, theta)[()]

    def db_of_psi_dtheta(self, psi: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $dB(\psi, \theta)/d\theta$, in Normalized Units."""
        return self._db_of_psi_dtheta(psi, theta)[()]

    def db_of_psip_dtheta(self, psip: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $dB(\psi_p, \theta)/d\theta$, in Normalized Units."""
        return self._db_of_psip_dtheta(psip, theta)[()]


class Geometry(_ReprStrImpl):
    """Geometry related evaluation methods."""

    _r: _PyGeometry

    def __init__(self) -> None:
        self._r_of_psi = np.vectorize(self._r.r_of_psi)
        self._r_of_psip = np.vectorize(self._r.r_of_psip)
        self._psi_of_r = np.vectorize(self._r.psi_of_r)
        self._psip_of_r = np.vectorize(self._r.psip_of_r)
        self._rlab_of_psi = np.vectorize(self._r.rlab_of_psi)
        self._rlab_of_psip = np.vectorize(self._r.rlab_of_psip)
        self._zlab_of_psi = np.vectorize(self._r.rlab_of_psi)
        self._zlab_of_psip = np.vectorize(self._r.rlab_of_psip)
        self._jacobian_of_psi = np.vectorize(self._r.rlab_of_psi)
        self._jacobian_of_psip = np.vectorize(self._r.rlab_of_psip)

    @property
    def baxis(self) -> float:
        r"""The magnetic field strength on the axis $B_0$ in $[T]$."""
        return self._r.baxis

    @property
    def raxis(self) -> float:
        r"""The horizontal position of the magnetic axis $R_0$ in $[m]$."""
        return self._r.raxis

    @property
    def zaxis(self) -> float:
        r"""The vertical position of the magnetic axis in $[m]$."""
        return self._r.zaxis

    @property
    def rgeo(self) -> float:
        r"""The horizontal position of the geometric axis (device major radius) in $[m]$."""
        return self._r.rgeo

    @property
    def rlast(self) -> float:
        r"""The $r$ coordinate's value at the last closed flux surface in $[m]$."""
        return self._r.rgeo

    @property
    def psi_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{LCFS}$."""
        return self._r.psi_last

    @property
    def psip_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{p,LCFS}$."""
        return self._r.psip_last

    @property
    def rlab_last(self) -> Array1:
        r"""The last $R$ values that correspond to the device's last closed flux surface, in $[m]$."""
        return self._r.rlab_last

    @property
    def zlab_last(self) -> Array1:
        r"""The last $Z$ values that correspond to the device's last closed flux surface, in $[m]$."""
        return self._r.zlab_last

    def r_of_psi(self, psi: ArrayLike) -> Array:
        r"""Calculates $r(\psi)$, where $\psi$ in Normalized Units and $r$ in $[m]$."""
        return self._r_of_psi(psi)[()]

    def r_of_psip(self, psip: ArrayLike) -> Array:
        r"""Calculates $r(\psi_p)$, where $\psi_p$ in Normalized Units and $r$ in $[m]$."""
        return self._r_of_psip(psip)[()]

    def psi_of_r(self, r: ArrayLike) -> Array:
        r"""Calculates $\psi(r)$, where $\psi$ in Normalized Units and $r$ in $[m]$."""
        return self._psi_of_r(r)[()]

    def psip_of_r(self, r: ArrayLike) -> Array:
        r"""Calculates $\psi_p(r)$, where $\psi_p$ in Normalized Units and $r$ in $[m]$."""
        return self._psip_of_r(r)[()]

    def rlab_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $R(\psi, \theta)$, where $\psi$ in Normalized Units and $R$ in $[m]$."""
        return self._rlab_of_psi(psi, theta)[()]

    def rlab_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $R(\psi_p, \theta)$, where $\psi_p$ in Normalized Units and $R$ in $[m]$."""
        return self._rlab_of_psip(psip, theta)[()]

    def zlab_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $Z(\psi, \theta)$, where $\psi$ in Normalized Units and $R$ in $[m]$."""
        return self._zlab_of_psi(psi, theta)[()]

    def zlab_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates $Z(\psi_p, \theta)$, where $\psi_p$ in Normalized Units and $R$ in $[m]$."""
        return self._zlab_of_psip(psip, theta)[()]

    def jacobian_of_psi(self, psi: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates the Jacobian $J(\psi, \theta)$, where $\psi$ in Normalized Units and $R$ in $[m]$."""
        return self._jacobian_of_psi(psi, theta)[()]

    def jacobian_of_psip(self, psip: ArrayLike, theta: ArrayLike) -> Array:
        r"""Calculates the Jacobian $R(\psi_p, \theta)$, where $\psi_p$ in Normalized Units and $R$ in $[m]$."""
        return self._jacobian_of_psip(psip, theta)[()]


class Mode(_ReprStrImpl):
    r"""Single perturbation mode related evaluation methods."""

    _r: _PyMode

    def __init__(self) -> None:
        self._ampl_of_psi = np.vectorize(self._r.ampl_of_psi)
        self._ampl_of_psip = np.vectorize(self._r.ampl_of_psip)
        self._phase_of_psi = np.vectorize(self._r.phase_of_psi)
        self._phase_of_psip = np.vectorize(self._r.phase_of_psip)
        self._m_of_psi = np.vectorize(self._r.m_of_psi)
        self._m_of_psip = np.vectorize(self._r.m_of_psi)
        self._dm_dpsi = np.vectorize(self._r.dm_dpsi)
        self._dm_dpsip = np.vectorize(self._r.dm_dpsip)
        self._dm_of_psi_dtheta = np.vectorize(self._r.dm_of_psi_dtheta)
        self._dm_of_psip_dtheta = np.vectorize(self._r.dm_of_psip_dtheta)
        self._dm_of_psi_dzeta = np.vectorize(self._r.dm_of_psi_dzeta)
        self._dm_of_psip_dzeta = np.vectorize(self._r.dm_of_psip_dzeta)
        self._dm_of_psi_dt = np.vectorize(self._r.dm_of_psi_dt)
        self._dm_of_psip_dt = np.vectorize(self._r.dm_of_psip_dt)

    @property
    def psi_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{LCFS}$."""
        return self._r.psi_last

    @property
    def psip_last(self) -> float:
        r"""The value of the last closed toroidal flux $\psi_{p,LCFS}$."""
        return self._r.psip_last

    @property
    def m(self) -> int:
        r"""The poloidal mode number $m$."""
        return self._r.m

    @property
    def n(self) -> int:
        r"""The toroidal mode number $n$."""
        return self._r.n

    def ampl_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates **the amplitude** $\alpha(\psi, \theta, \zeta, t)$, in Normalized Units."""
        return self._ampl_of_psi(psi, theta, zeta, t)[()]

    def ampl_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates **the amplitude** $\alpha(\psi_p, \theta, \zeta, t)$, in Normalized Units."""
        return self._ampl_of_psip(psip, theta, zeta, t)[()]

    def phase_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates $\phi(\psi, \theta, \zeta, t)$, in Normalized Units."""
        return self._phase_of_psi(psi, theta, zeta, t)[()]

    def phase_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates $\phi(\psi_p, \theta, \zeta, t)$, in Normalized Units."""
        return self._phase_of_psip(psip, theta, zeta, t)[()]

    def m_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's value $m(\psi, \theta, \zeta, t)$, in Normalized Units."""
        return self._m_of_psi(psi, theta, zeta, t)[()]

    def m_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's value $m(\psi_p, \theta, \zeta, t)$, in Normalized Units."""
        return self._m_of_psip(psip, theta, zeta, t)[()]

    def dm_dpsi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/d\psi$, in Normalized Units."""
        return self._dm_dpsi(psi, theta, zeta, t)[()]

    def dm_dpsip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi_p, \theta, \zeta, t)/d\psi_p$, in Normalized Units."""
        return self._dm_dpsip(psip, theta, zeta, t)[()]

    def dm_of_psi_dtheta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/d\theta$, in Normalized Units."""
        return self._dm_of_psi_dtheta(psi, theta, zeta, t)[()]

    def dm_of_psip_dtheta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi_p, \theta, \zeta, t)/d\theta$, in Normalized Units."""
        return self._dm_of_psip_dtheta(psip, theta, zeta, t)[()]

    def dm_of_psi_dzeta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/d\zeta$, in Normalized Units."""
        return self._dm_of_psi_dzeta(psi, theta, zeta, t)[()]

    def dm_of_psip_dzeta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi_p, \theta, \zeta, t)/d\zeta$, in Normalized Units."""
        return self._dm_of_psip_dzeta(psip, theta, zeta, t)[()]

    def dm_of_psi_dt(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi, \theta, \zeta, t)/dt$, in Normalized Units."""
        return self._dm_of_psi_dt(psi, theta, zeta, t)[()]

    def dm_of_psip_dt(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> Array:
        r"""Calculates the mode's derivative $dm(\psi_p, \theta, \zeta, t)/dt$, in Normalized Units."""
        return self._dm_of_psip_dt(psip, theta, zeta, t)[()]
