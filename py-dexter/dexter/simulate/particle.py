from collections.abc import Sequence

from dexter.machine.machine import Machine
from dexter.simulate.initial import InitialConditions

from dexter.types import (
    Array1,
    EnergyPzetaPosition,
    IntegrationStatus,
    Intersection,
    OrbitType,
    SteppingMethod,
)

from dexter._core import _PyParticle, _PySolverParams, _PyIntersectParams
from dexter._utils import _ReprStrImpl


class Particle(_ReprStrImpl):
    r"""A Particle.

    By taking $\mu = 0$ and $\rho \rightarrow 0$, the particle traces magnetic field
    lines.

    Parameters
    ----------
    initial_conditions
        The [InitialConditions][dexter.InitialConditions] set.

    Example
    -------
    ```python title="Particle creation from a Boozer Coordinates set"
    >>> initial_conditions = dex.InitialConditions.Boozer(
    ...     t0=0,
    ...     flux0=dex.InitialFlux.Toroidal(0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     rho0=1e-4,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```
    ```python title="Particle creation from a Mixed Coordinates set"
    >>> initial_conditions = dex.InitialConditions.Mixed(
    ...     t0=0,
    ...     flux0=dex.InitialFlux.Poloidal(0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     pzeta0=-0.025,
    ...     mu0=7e-6,
    ... )
    >>> particle = dex.Particle(initial_conditions)

    ```
    """

    _r: _PyParticle

    def __init__(self, initial: InitialConditions) -> None:
        self._r = _PyParticle(initial._r)
        pass

    def integrate(
        self,
        machine: Machine,
        teval: tuple[float, float],
        *,
        stepping_method: SteppingMethod | None = "EnergyAdaptiveStep",
        max_steps: int | None = 1_000_000,
        first_step: float | None = 1e-1,
        safety_factor: float | None = 0.9,
        energy_rel_tol: float | None = 1e-12,
        energy_abs_tol: float | None = 1e-14,
        error_rel_tol: float | None = 1e-12,
        error_abs_tol: float | None = 1e-14,
    ):
        r"""Integrates the particle for a specific time interval.

        The time interval is in Normalized Units (inverse gyro-frequency).

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.
        teval
            The time span $(t_0, t_f)$ in which to integrate the particle, in Normalized Units.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Particle integration"
        >>> # Machine setup
        >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation=dex.Perturbation(
        ...         [
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=1, n=3, phase=0),
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>>
        >>> # Initial conditions setup
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux.Toroidal(0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>>
        >>> # Particle setup and integration
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.integrate(
        ...     machine=machine,
        ...     teval=(0, 1e2),
        ...     energy_rel_tol=1e-11,
        ...     energy_abs_tol=1e-13,
        ... )

        ```
        """
        solver_params = _PySolverParams(
            stepping_method=stepping_method,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
            energy_rel_tol=energy_rel_tol,
            energy_abs_tol=energy_abs_tol,
            error_rel_tol=error_rel_tol,
            error_abs_tol=error_abs_tol,
        )
        self._r.integrate(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            teval=teval,
            solver_params=solver_params,
        )

    def intersect(
        self,
        machine: Machine,
        intersection: Intersection,
        angle: float,
        turns: int,
        *,
        stepping_method: SteppingMethod | None = "EnergyAdaptiveStep",
        max_steps: int | None = 1_000_000,
        first_step: float | None = 1e-1,
        safety_factor: float | None = 0.9,
        energy_rel_tol: float | None = 1e-12,
        energy_abs_tol: float | None = 1e-14,
        error_rel_tol: float | None = 1e-12,
        error_abs_tol: float | None = 1e-14,
    ):
        r"""Integrates the particle, calculating its intersections with a constant $\theta$ or $\zeta$ surface.

        Using the method described by Hénon we can force the solver to step exactly on the intersection surface.

        The differences between two consecutive values of the corresponding angle variable are guaranteed to
        be $2\pi \pm \epsilon$, where $\epsilon$ a number smaller than the solver’s relative tolerance.

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.
        intersection
            The surface of section Σ, defined by an equation $\chi_i = \alpha$, where $\chi_i = \theta$ or
            $\zeta$.
        angle
            The constant that defines the surface of section.
        turns
            The number of intersections to calculate.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Particle intersection integration"
        >>> # Machine setup
        >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ...     perturbation=dex.Perturbation(
        ...         [
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=1, n=3, phase=0),
        ...             dex.FluteMode(epsilon=1e-3, lcfs=LCFS, m=2, n=3, phase=0),
        ...         ]
        ...     )
        ... )
        >>>
        >>> # Initial conditions and Intersection Parameters setup
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux.Toroidal(0.3),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-3,
        ...     mu0=7e-5,
        ... )
        >>>
        >>> # Particle setup and intersection
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.intersect(
        ...     machine=machine,
        ...     intersection="ConstZeta",
        ...     angle=3.1415,
        ...     turns=5,
        ... )

        ```
        """
        solver_params = _PySolverParams(
            stepping_method=stepping_method,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
            energy_rel_tol=energy_rel_tol,
            energy_abs_tol=energy_abs_tol,
            error_rel_tol=error_rel_tol,
            error_abs_tol=error_abs_tol,
        )
        intersect_params = _PyIntersectParams(
            intersection=intersection,
            angle=angle,
            turns=turns,
        )

        self._r.intersect(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            intersect_params=intersect_params,
            solver_params=solver_params,
        )

    def close(
        self,
        machine: Machine,
        periods: int | None = 1,
        *,
        stepping_method: SteppingMethod | None = "EnergyAdaptiveStep",
        max_steps: int | None = 1_000_000,
        first_step: float | None = 1e-1,
        safety_factor: float | None = 0.9,
        energy_rel_tol: float | None = 1e-12,
        energy_abs_tol: float | None = 1e-14,
        error_rel_tol: float | None = 1e-12,
        error_abs_tol: float | None = 1e-14,
    ):
        r"""Integrates the particle for a certain amount of $\theta-\psi$ periods.

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.
        periods
            The amount of periods to integrate.

        Other Parameters
        ----------------
        stepping_method
            The optimal step calculation method. Defaults to "EnergyAdaptiveStep".
        max_steps
            The maximum amount of steps a particle can make before terminating its integration. Defaults to
            1.000.000.
        first_step
            The initial time step for the RKF45 adaptive step method. The value is empirical. Defaults to
            1e-1.
        safety_factor
            The safety factor of the solver. Should be less than 1.0. Defaults to 0.9.
        energy_rel_tol
            The relative tolerance of the energy difference in every step. Defaults to 1e-12.
        energy_abs_tol
            The absolute tolerance of the energy difference in every step. Defaults to 1e-14.
        error_rel_tol
            The relative tolerance of the local truncation error in every step. Defaults to 1e-12.
        error_abs_tol
            The absolute tolerance of the local truncation error in every step. Defaults to 1e-14.

        Example
        -------
        ```python title="Particle orbit closing"
        >>> # Machine setup
        >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> # Initial conditions setup
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux.Toroidal(0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-5,
        ...     mu0=7e-6,
        ... )
        >>>
        >>> # Particle setup and integration
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.close(machine)

        ```
        """
        solver_params = _PySolverParams(
            stepping_method=stepping_method,
            max_steps=max_steps,
            first_step=first_step,
            safety_factor=safety_factor,
            energy_rel_tol=energy_rel_tol,
            energy_abs_tol=energy_abs_tol,
            error_rel_tol=error_rel_tol,
            error_abs_tol=error_abs_tol,
        )

        self._r.close(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
            perturbation=machine.perturbation._r,
            periods=periods if periods is not None else 1,
            solver_params=solver_params,
        )

    def classify(
        self,
        machine: Machine,
    ):
        r"""Classifies the particle’s orbit using its position on the $(E, P_\zeta, \mu=const)$
        plane without integrating.

        This routine only works for LAR-like equilibria.

        Parameters
        ----------
        machine
            The machine in which to integrate the particle.

        Example
        -------
        ```python title="Particle classification"
        >>> # Machine setup
        >>> LCFS = dex.LastClosedFluxSurface.Toroidal(0.45)
        >>> machine = dex.Machine(
        ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.1, lcfs=LCFS),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> # Initial conditions setup
        >>> initial_conditions = dex.InitialConditions.Boozer(
        ...     t0=0,
        ...     flux0=dex.InitialFlux.Toroidal(0.1),
        ...     theta0=3.14,
        ...     zeta0=0,
        ...     rho0=1e-4,
        ...     mu0=7e-6,
        ... )
        >>>
        >>> # Particle setup and classification
        >>> particle = dex.Particle(initial_conditions)
        >>> particle.classify(machine)

        ```
        """
        self._r.classify(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
        )

    @property
    def initial_conditions(self) -> InitialConditions:
        """The Particle's [`InitialConditions`][dexter.InitialConditions]."""
        return InitialConditions._wrap(self._r.initial_conditions)

    @property
    def integration_status(self) -> IntegrationStatus:
        """The Particle's [`IntegrationStatus`][dexter.IntegrationStatus]."""
        return self._r.integration_status

    @property
    def steps_taken(self) -> int:
        """The total number of steps taken during the integration.

        This number is not necessarily the same as the number of steps stored."""
        return self._r.steps_taken

    @property
    def steps_stored(self) -> int:
        """The total number of steps stored in the time arrays."""
        return self._r.steps_stored

    @property
    def duration(self) -> str:
        """The duration of the integration routine in Normalized Units."""
        return self._r.duration

    @property
    def initial_energy(self) -> float:
        """The Particle's initial energy in Normalized Units."""
        if self._r.initial_energy is None:
            raise AttributeError("`initial_energy` has not been calculated")
        return self._r.initial_energy

    @property
    def final_energy(self) -> float:
        """The Particle's final energy in Normalized Units."""
        if self._r.final_energy is None:
            raise AttributeError("`final_energy` has not been calculated")
        return self._r.final_energy

    @property
    def energy_var(self) -> float:
        """The variance of the Particle's energy array in Normalized Units."""
        if self._r.energy_var is None:
            raise AttributeError("`energy_var` has not been calculated")
        return self._r.energy_var

    @property
    def energy_pzeta_position(self) -> EnergyPzetaPosition:
        """The Particle's [`EnergyPzetaPosition`][dexter.EnergyPzetaPosition]."""
        if self._r.energy_pzeta_position is None:
            raise AttributeError("`energy_pzeta_position` has not been calculated")
        return self._r.energy_pzeta_position

    @property
    def orbit_type(self) -> OrbitType:
        """The Particle's [`OrbitType`][dexter.OrbitType]."""
        if self._r.orbit_type is None:
            raise AttributeError("`orbit_rype` has not been calculated")
        return self._r.orbit_type

    @property
    def omega_theta(self) -> float:
        r"""The Particle's $\omega_\theta$ frequency in Normalized Units."""
        if self._r.omega_theta is None:
            raise AttributeError("`omega_rheta` has not been calculated")
        return self._r.omega_theta

    @property
    def omega_zeta(self) -> float:
        r"""The Particle's $\omega_\zeta$ frequency in Normalized Units."""
        if self._r.omega_zeta is None:
            raise AttributeError("`omega_zeta` has not been calculated")
        return self._r.omega_zeta

    @property
    def qkinetic(self) -> float:
        r"""The Particle's $q_{kin}$."""
        if self._r.qkinetic is None:
            raise AttributeError("`qkinetic` has not been calculated")
        return self._r.qkinetic

    def print_caches(self):
        r"""Prints the Particle's accelerators' hits and misses."""
        self._r.print_caches()

    def discard_arrays(self):
        """Discards the time series arrays."""
        self._r.discard_arrays()

    @property
    def flux_cache_hits(self) -> int:
        """The magnetic flux' Accelerator cache hits."""
        return self._r.flux_cache_hits

    @property
    def flux_cache_misses(self) -> int:
        """The magnetic flux' Accelerator cache misses."""
        return self._r.flux_cache_misses

    @property
    def theta_cache_hits(self) -> int:
        r"""The $\theta$ coordinate's Accelerator cache misses."""
        return self._r.theta_cache_hits

    @property
    def theta_cache_misses(self) -> int:
        r"""The $\theta$ coordinate's Accelerator cache misses."""
        return self._r.theta_cache_misses

    @property
    def mode_cache_hits(self) -> int:
        """The modes' Accelerator cache hits."""
        return self._r.mode_cache_hits

    @property
    def mode_cache_misses(self) -> int:
        """The modes' Accelerator cache misses."""
        return self._r.mode_cache_misses

    @property
    def t_array(self) -> Array1:
        r"""The $t$ array."""
        return self._r.get_array("t_array")

    @property
    def psi_array(self) -> Array1:
        r"""The $\psi$ array."""
        return self._r.get_array("psi_array")

    @property
    def psip_array(self) -> Array1:
        r"""The $\psi_p$ array."""
        return self._r.get_array("psip_array")

    @property
    def theta_array(self) -> Array1:
        r"""The $\theta$ array."""
        return self._r.get_array("theta_array")

    @property
    def zeta_array(self) -> Array1:
        r"""The $\zeta$ array."""
        return self._r.get_array("zeta_array")

    @property
    def rho_array(self) -> Array1:
        r"""The $\rho_{||}$ array."""
        return self._r.get_array("rho_array")

    @property
    def mu_array(self) -> Array1:
        r"""The $\mu$ array."""
        return self._r.get_array("mu_array")

    @property
    def ptheta_array(self) -> Array1:
        r"""The $P_\theta$ array."""
        return self._r.get_array("ptheta_array")

    @property
    def pzeta_array(self) -> Array1:
        r"""The $P_\zeta$ array."""
        return self._r.get_array("pzeta_array")

    @property
    def energy_array(self) -> Array1:
        r"""The energy array."""
        return self._r.get_array("energy_array")
