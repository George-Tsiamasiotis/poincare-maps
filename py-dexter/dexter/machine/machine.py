"""Defines `Machine`, a container for all information a device."""

from dexter.types import ArrayLike, ParticleSpecies, Unit
from pint.facets.plain import PlainQuantity

from dexter.machine.bfields import BfieldObject
from dexter.machine.currents import CurrentObject
from dexter.machine.geometries import GeometryObject
from dexter.machine.perturbation import Perturbation
from dexter.machine.qfactors import QfactorObject

from dexter.machine._registry import _Registry


class Machine:
    r"""A Machine.

    Container for all information about the magnetic field, geometry and perturbations of a device.

    Parameters
    ----------
    qfactor
        The machine's's q-factor.
    current
        The machine's's current.
    bfield
        The machine's's bfield.
    perturbation
        The machine's's perturbation.
    geometry
        The machine's's geometry. This field is not required for routines, but is required for
        conversion to laboratory coordinates and SI units.
    species
        The particle species under study. This only affects unit conversions to SI.

    Example
    -------

    ```python title="Machine creation"
    >>> lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    >>> machine = dex.Machine(
    ...     qfactor=dex.UnityQfactor(lcfs),
    ...     current=dex.LarCurrent(),
    ...     bfield=dex.LarBfield(),
    ...     perturbation=dex.Perturbation(
    ...          [
    ...             dex.FluteMode(1e-5, lcfs, 2, 1, 0),
    ...             dex.FluteMode(1e-5, lcfs, 3, 1, 0),
    ...          ]
    ...     ),
    ... )

    ```

    """

    _geometry: GeometryObject | None
    _qfactor: QfactorObject
    _current: CurrentObject
    _bfield: BfieldObject
    _perturbation: Perturbation
    _reg = _Registry

    def __init__(
        self,
        qfactor: QfactorObject,
        current: CurrentObject,
        bfield: BfieldObject,
        *,
        geometry: GeometryObject | None = None,
        perturbation: Perturbation = Perturbation([]),
        species: ParticleSpecies = "Proton",
    ) -> None:
        self._qfactor = qfactor
        self._current = current
        self._bfield = bfield
        self._geometry = geometry
        self._perturbation = perturbation

        if geometry is not None:
            self._reg = _Registry(case_sensitive=False, cache_folder=":auto:")
            self._reg.define_normalizations(geometry.raxis, geometry.baxis, species)
        else:
            self._reg = None

    def quantity(self, value: float | ArrayLike, units: Unit) -> PlainQuantity:
        if self._reg is None:
            raise RuntimeError("UnitRegistry has not been defined")
        else:
            return self._reg.Quantity(value, units)

    @property
    def geometry(self) -> GeometryObject:
        """The machine's ['GeometryObject'](dexter.GeometryObject)."""
        if self._geometry is None:
            raise AttributeError("`GeometryObject` has not been defined")
        else:
            return self._geometry

    @property
    def qfactor(self) -> QfactorObject:
        """The machine's ['QfactorObject'](dexter.QfactorObject)."""
        return self._qfactor

    @property
    def current(self) -> CurrentObject:
        """The machine's ['CurrentObject'](dexter.CurrentObject)."""
        return self._current

    @property
    def bfield(self) -> BfieldObject:
        """The machine's ['BfieldObject'](dexter.BfieldObject)."""
        return self._bfield

    @property
    def perturbation(self) -> Perturbation:
        """The machine's 'Perturbation'."""
        return self._perturbation

    @perturbation.setter
    def perturbation(self, perturbation: Perturbation):
        """Sets the machine's's Perturbation."""
        self._perturbation = perturbation

    @property
    def psi_last(self) -> float:
        r"""The value of the last closed toroidal flux surface, $\psi_{LCFS}$."""
        return self.qfactor.psi_last

    @property
    def psip_last(self) -> float:
        r"""The value of the last closed poloidal flux surface, $\psi_{p,LCFS}$."""
        return self.qfactor.psip_last

    @property
    def baxis(self) -> float:
        r"""The magnetic field strength on the axis, $B_{axis}$."""
        if self._geometry is not None:
            return self._geometry.baxis
        else:
            raise AttributeError("`baxis` is not defined")

    @property
    def raxis(self) -> float:
        r"""The device's major radius, $R_{axis}$."""
        if self._geometry is not None:
            return self._geometry.raxis
        else:
            raise AttributeError("`raxis` is not defined")

    @property
    def zaxis(self) -> float:
        r"""The vertical position of the magnetic axis in $[m]$."""
        if self._geometry is not None:
            return self._geometry.zaxis
        else:
            raise AttributeError("`zaxis` is not defined")

    @property
    def rgeo(self) -> float:
        r"""The horizontal position of the geometric axis (device major radius) in $[m]$."""
        if self._geometry is not None:
            return self._geometry.rgeo
        else:
            raise AttributeError("`rgeo` is not defined")

    @property
    def rlast(self) -> float:
        r"""The radial coordinate's value at the last closed flux surface, $r_{LCFS}$."""
        if self._geometry is not None:
            return self._geometry.rlast
        else:
            raise AttributeError("`rlast` is not defined")

    def __str__(self) -> str:
        string = "Machine:\n"
        if self._geometry is not None:
            string += str(getattr(self, "geometry", "")) + "\n"
        string += str(getattr(self, "qfactor", "")) + "\n"
        string += str(getattr(self, "current", "")) + "\n"
        string += str(getattr(self, "bfield", "")) + "\n"
        string += str(getattr(self, "perturbation", "")) + "\n"
        return string

    def __repr__(self) -> str:
        return self.__str__()
