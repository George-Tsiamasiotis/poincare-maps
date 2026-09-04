"""Defines types associated with Particle and Queue initialization."""

from dexter.types import FluxCoordinate, CoordinateSet

from dexter._utils import _ReprStrImpl
from dexter._core import _PyInitialFlux, _PyInitialConditions


class InitialFlux(_ReprStrImpl):
    """Defines the flux coordinate initial value.

    This type is instantiated through the [`Toroidal`][dexter.InitialFlux.Toroidal] and
    [`Poloidal`][dexter.InitialFlux.Poloidal] class methods.
    """

    _r: _PyInitialFlux

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def Toroidal(cls, value: float) -> InitialFlux:
        r"""Defines the initial value with respect to $\psi$.

        Parameters
        ----------
        value
            The value of the toroidal magnetic flux at the last closed flux surface.

        Example
        -------
        ```python title="InitialFlux creation"
        >>> flux0 = dex.InitialFlux.Toroidal(0.05)

        ```
        """
        obj = InitialFlux.__new__(InitialFlux)
        obj._r = _PyInitialFlux.toroidal(value)
        return obj

    @classmethod
    def Poloidal(cls, value: float) -> InitialFlux:
        r"""Defines the initial value with respect to $\psi_p$.

        Parameters
        ----------
        value
            The value of the poloidal magnetic flux at the last closed flux surface.

        Example
        -------
        ```python title="InitialFlux creation"
        >>> flux0 = dex.InitialFlux.Poloidal(0.05)

        ```
        """
        obj = InitialFlux.__new__(InitialFlux)
        obj._r = _PyInitialFlux.poloidal(value)
        return obj

    @property
    def value(self) -> float:
        """The value of the magnetic flux."""
        return self._r.value

    @property
    def kind(self) -> FluxCoordinate:
        """The kind of the magnetic flux."""
        return self._r.kind


class InitialConditions(_ReprStrImpl):
    r"""Initial conditions set for a Particle.

    This type is instantiated through the [`Boozer`][dexter.InitialConditions.Boozer] and
    [`Mixed`][dexter.InitialConditions.Mixed] class methods.
    """

    _r: _PyInitialConditions

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def Boozer(
        cls,
        t0: float,
        flux0: InitialFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ) -> InitialConditions:
        r"""Creates initial conditions for a Particle in Boozer coordinates.

        The initial conditions are defined on the
        $(t, \psi, \theta, \zeta, \rho, \mu)$ or
        $(t, \psi_p, \theta, \zeta, \rho, \mu)$
        space, depending on the value of `flux0`.

        Parameters
        ----------
        t0
            The initial time, in Normalized Units.
        flux0
            The initial $\psi / \psi_p$, in Normalized Units.
        theta0
            The initial $\theta$ angle, in rads.
        zeta0
            The initial $\zeta$ angle, in rads.
        rho0
            The initial $\rho_{||}$, in Normalized Units.
        mu0
            The initial magnetic moment $\mu$, in Normalized Units.

        Example
        -------
        ```python title="InitialConditions definition in Boozer-Toroidal coordinates"
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux.Toroidal(0.01),  # ψ0 = 0.01
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )

        ```
        """
        obj = InitialConditions.__new__(InitialConditions)
        obj._r = _PyInitialConditions.boozer(t0, flux0._r, theta0, zeta0, rho0, mu0)
        return obj

    @classmethod
    def Mixed(
        cls,
        t0: float,
        flux0: InitialFlux,
        theta0: float,
        zeta0: float,
        pzeta0: float,
        mu0: float,
    ) -> InitialConditions:
        r"""Creates initial conditions for a Particle in Mixed coordinates.

        The initial conditions are defined on the
        $(t, \psi, \theta, \zeta, P_\zeta, \mu)$ or
        $(t, \psi_p, \theta, \zeta, P_\zeta, \mu)$
        space, depending on the value of `flux0`.

        Parameters
        ----------
        t0
            The initial time, in Normalized Units.
        flux0
            The initial $\psi / \psi_p$, in Normalized Units.
        theta0
            The initial $\theta$ angle, in rads.
        zeta0
            The initial $\zeta$ angle, in rads.
        pzeta0
            The initial $P_\zeta$, in Normalized Units.
        mu0
            The initial magnetic moment $\mu$, in Normalized Units.

        Example
        -------
        ```python title="InitialConditions definition in Mixed-Poloidal coordinates"
        >>> flux0=dex.InitialFlux.Poloidal(0.02)  # ψ0 = 0.02
        >>> initial_conditions = dex.InitialConditions.Mixed(
        ...     t0=0,
        ...     flux0=flux0,
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     pzeta0=-0.4 * flux0.value,
        ...     mu0=7e-6,
        ... )

        ```
        """
        obj = InitialConditions.__new__(InitialConditions)
        obj._r = _PyInitialConditions.mixed(t0, flux0._r, theta0, zeta0, pzeta0, mu0)
        return obj

    @property
    def t0(self) -> float:
        """The initial time, in Normalized Units."""
        return self._r.t0

    @property
    def flux0(self) -> InitialFlux:
        r"""The initial $\psi / \psi_p$, in Normalized Units."""
        initial_flux = InitialFlux.__new__(InitialFlux)
        initial_flux._r = self._r.flux0
        return initial_flux

    @property
    def theta0(self) -> float:
        r"""The initial $\theta$ angle, in Normalized Units."""
        return self._r.theta0

    @property
    def zeta0(self) -> float:
        r"""The initial $\zeta$ angle, in Normalized Units."""
        return self._r.zeta0

    @property
    def rho0(self) -> float:
        r"""The initial parallel radius $\rho_{||}$, in Normalized Units."""
        if self._r.rho0 is None:
            raise AttributeError("'rho0' has not been defined")
        else:
            return self._r.rho0

    @property
    def pzeta0(self) -> float:
        r"""The initial canonical momentum $P_\zeta$, in Normalized Units."""
        if self._r.pzeta0 is None:
            raise AttributeError("'pzeta0' has not been defined")
        else:
            return self._r.pzeta0

    @property
    def mu0(self) -> float:
        r"""The initial magnetic moment $\mu$, in Normalized Units."""
        return self._r.mu0

    @property
    def coordinate_set(self) -> CoordinateSet:
        """The kind of InitialConditions set."""
        return self._r.coordinate_set
