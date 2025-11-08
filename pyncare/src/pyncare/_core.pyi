import numpy as np

class Qfactor:
    """q-factor reconstructed from a NetCDF file.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psi_wall: float
        The value of the toroidal flux ψ at the wall.
    psip_data: np.ndarray
        The NetCDF ψp data used to construct the q(ψp) and ψ(ψp) splines.
    q_data: np.ndarray
        The NetCDF q data used to construct the q(ψp) spline.
    psi_data: np.ndarray
        The NetCDF ψp data used to construct the ψ(ψp) spline.
    q_data_derived: np.ndarray
        The q values, as calculated from dψ/dψp, at the ψp data.
    """

    path: str
    typ: str
    psip_wall: float
    psi_wall: float
    psip_data: np.ndarray
    q_data: np.ndarray
    psi_data: np.ndarray
    q_data_derived: np.ndarray

    def __init__(self, path: str, typ: str) -> None:
        """q-factor reconstructed from a NetCDF file.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: str
            The type of Interpolation. Available types: "Linear", "Cubic", "Akima", "AkimaPeriodic",
            "Steffen".
        """

    def q(self, psip: float) -> float:
        """The q value evaluated at ψp"""

    def psi(self, psip: float) -> float:
        """The ψ value evaluated at ψp"""

class Currents:
    """Plasma current reconstructed from a NetCDF file.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psip_data: np.ndarray
        The NetCDF ψp data used to construct the g(ψp) and I(ψp) splines.
    g_data: np.ndarray
        The NetCDF g data used to construct the g(ψp) spline.
    i_data: np.ndarray
        The NetCDF I data used to construct the I(ψp) spline.
    """

    path: str
    typ: str
    psip_wall: float
    psip_data: np.ndarray
    g_data: np.ndarray
    i_data: np.ndarray

    def __init__(self, path: str, typ: str):
        """Plasma current reconstructed from a NetCDF file.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: str
            The type of Interpolation. Available types: "Linear", "Cubic", "Akima", "AkimaPeriodic",
            "Steffen".
        """

    def g(self, psip: float) -> float:
        """The g value evaluated at ψp"""

    def i(self, psip: float) -> float:
        """The I value evaluated at ψp"""

    def dg_dpsip(self, psip: float) -> float:
        """The dg/dψp value evaluated at ψp"""

    def di_dpsip(self, psip: float) -> float:
        """The dI/dψp value evaluated at ψp"""

class Bfield:
    """Magnetic field reconstructed from a NetCDF file.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    baxis: float
        The magnetic field strength on the axis in [T].
    raxis: float
        The major radius in [m].
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psip_data: np.ndarray
        The NetCDF ψp data used to construct the b(ψp, θ) spline.
    theta_data: np.ndarray
        The NetCDF θ data used to construct the b(ψp, θ) spline.
    b_data: np.ndarray
        The NetCDF b data used to construct the b(ψp, θ) spline.
    r_data: np.ndarray
        The NetCDF R data used to construct the R(ψp, θ) spline.
    z_data: np.ndarray
        The NetCDF Z data used to construct the Z(ψp, θ) spline.
    db_dpsip_data: np.ndarray
        The db/dψp values evaluated at the (ψp, θ) data, through interpolation.
    db_dtheta_data: np.ndarray
        The db/dψp values evaluated at the (ψp, θ) data, through interpolation.
    """

    path: str
    typ: str
    baxis: float
    raxis: float
    psip_wall: float
    psip_data: np.ndarray
    theta_data: np.ndarray
    b_data: np.ndarray
    r_data: np.ndarray
    z_data: np.ndarray
    db_dpsip_data: np.ndarray
    db_dtheta_data: np.ndarray

    def __init__(self, path: str, typ: str):
        """Magnetic field reconstructed from a NetCDF file.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: str
            The type of Interpolation. Available types: "Bilinear", "Bicubic".
        """

    def b(self, psip: float, theta: float) -> float:
        """The b value evaluated at (ψp, θ)"""

    def db_dpsip(self, psip: float, theta: float) -> float:
        """The db/dψp value evaluated at (ψp, θ)"""

    def db_dtheta(self, psip: float, theta: float) -> float:
        """The db/dθ value evaluated at (ψp, θ)"""

    def d2b_dpsip2(self, psip: float, theta: float) -> float:
        """The d2b/dψp2 value evaluated at (ψp, θ)"""

    def d2b_dtheta2(self, psip: float, theta: float) -> float:
        """The d2b/dθ2 value evaluated at (ψp, θ)"""

    def d2b_dpsip_dtheta(self, psip: float, theta: float) -> float:
        """The d2b/dψpdθ value evaluated at (ψp, θ)"""

class Harmonic:
    """A single perturbation harmonic.

    Attributes
    ----------
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psip_data: float
        The NetCDF ψp data used to construct the a(ψp) spline.
    a_data: float
        The NetCDF a data used to construct the a(ψp) spline.
    """

    path: str
    typ: str
    m: float
    n: float
    phase: float
    psip_wall: float
    psip_data: np.ndarray
    a_data: np.ndarray

    def __init__(self, path: str, typ: str, m: float, n: float, phase: float):
        """Creates a single perturbation harmonic.

        Parameters
        ----------
        path: str
            The path to the NetCDF file.
        typ: str
            The type of Interpolation. Available types: "Linear", "Cubic", "Akima", "AkimaPeriodic",
            "Steffen".
        m: float
            The `θ` frequency number.
        n: float
            The `ζ` frequency number.
        phase: float
            The initial phase of the harmonic.
        """

    def h(self, psip: float, theta: float, zeta: float) -> float:
        """The h value (value of the whole harmonic) evaluated at (ψp, θ, ζ)."""

    def dh_dpsip(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dψp value evaluated at (ψp, θ, ζ)."""

    def dh_dtheta(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dθ value evaluated at (ψp, θ, ζ)."""

    def dh_dzeta(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dζ value evaluated at (ψp, θ, ζ)."""

    def dh_dt(self, psip: float, theta: float, zeta: float) -> float:
        """The dh/dt value evaluated at (ψp, θ, ζ)."""

class Perturbation:
    """A sum of different perturbation harmonics."""

    harmonics: list[Harmonic]

    def __init__(self, harmonics: list[Harmonic]):
        """Creates a Perturbation.

        Parameters
        ----------
        harmonics: list[Harmonics]
            The list of harmonics that appear in the perturbation.
        """

    def __getitem__(self, n: int) -> Harmonic:
        """Returns the n-th harmonic"""

    def p(self, psip: float, theta: float, zeta: float) -> float:
        """The p value (value of the whole harmonic) evaluated at (ψp, θ, ζ)."""

    def dp_dpsip(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dψp value evaluated at (ψp, θ, ζ)."""

    def dp_dtheta(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dθ value evaluated at (ψp, θ, ζ)."""

    def dp_dzeta(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dζ value evaluated at (ψp, θ, ζ)."""

    def dp_dt(self, psip: float, theta: float, zeta: float) -> float:
        """The dp/dt value evaluated at (ψp, θ, ζ)."""

# =========================================================================================

class InitialConditions:
    """A set of initial conditions."""

    # TODO:
    # t0: float
    # theta0: float
    # psip0: float
    # rho0: float
    # zeta0: float
    # mu: float

    def __init__(
        self,
        time0: float,
        theta0: float,
        psip0: float,
        rho0: float,
        zeta0: float,
        mu: float,
    ):
        """Creates a set of initial conditions.

        Parameters
        ----------
        time0: float
            The initial time.
        theta0: float
            The initial `θ` angle.
        psip0: float
            The initial poloidal magnetic flux `ψp`.
        rho0: float
            The initial parallel gyro radius `ρ`.
        zeta0: float
            The initial `ζ` angle.
        mu: float
            The magnetic moment `μ`.
        """

class Particle:
    """A particle."""

    evolution: Evolution
    status: str

    def __init__(self, initial: InitialConditions):
        """Creates a Particle from an `InitialConditions` set.

        Parameters
        ----------
        initial: InitialConditions
            The initial conditions set.
        """

    def integrate(
        self,
        qfactor: Qfactor,
        bfield: Bfield,
        currents: Currents,
        perturbation: Perturbation,
        t_eval: list[float],
    ):
        """Integrates the particle, storing its evolution.

        Parameters
        ----------
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        t_eval: list[float]
            The time span [t0, tf] in which to integrate the particle.
        """

    def map(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
        params: MappingParameters,
    ):
        """Integrates the particle, storing its intersections with the Poincare
        surface defined by `MappingParameters`.

        Parameters
        ----------
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        params: MappingParameters
            The parameters of the Poincare mapping.
        """

class Evolution:
    """Time series of the particle's orbit.

    Not meant to be constructed. It is stored as a particle's attribute.
    """

    time: np.ndarray
    theta: np.ndarray
    psip: np.ndarray
    rho: np.ndarray
    zeta: np.ndarray
    psi: np.ndarray
    ptheta: np.ndarray
    pzeta: np.ndarray
    steps_taken: int
    steps_stored: int

class MappingParameters:
    """Defines all the necessary parameters of a Poincare Map.

    Attributes
    ----------
    section: str
        The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    alpha: float
        The constant that defines the surface of section (modulo 2π).
    intersections: int
        The number of interections to calculate.
    """

    section: str
    alpha: float
    intersections: int

    def __init__(self, section: str, alpha: float, intersections: int):
        """Defines all the necessary parameters of a Poincare Map.

        Parameters
        ----------
        section: str
            The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
        alpha: float
            The constant that defines the surface of section (modulo 2π).
        intersections: int
            The number of interections to calculate.
        """

# =========================================================================================

class PoincareInit:
    """Sets up the initial conditions for a poincare map."""

    thetas: np.ndarray
    psips: np.ndarray
    rhos: np.ndarray
    zetas: np.ndarray
    mus: np.ndarray

    def __init__(
        self,
        thetas: np.ndarray,
        psips: np.ndarray,
        rhos: np.ndarray,
        zetas: np.ndarray,
        mus: np.ndarray,
    ):
        """Constructor.

        Parameters
        ----------
        thetas: np.ndarray
            The initial `θ` angles.
        psips: np.ndarray
            The initial poloidal magnetic fluxes `ψp`.
        rhos: np.ndarray
            The initial parallel gyro radii `ρ`.
        zetas: np.ndarray
            The initial `ζ` angles.
        mus: np.ndarray
            The magnetic moments `μ`.
        """

class Poincare:
    """Calculates Poincare maps."""

    init: PoincareInit
    section: str
    alpha: int
    intersection: int
    angles: np.ndarray
    fluxes: np.ndarray

    def __init__(self, init: PoincareInit, params: MappingParameters):
        """Constructor

        Parameters
        ----------
        init: PoincareInit
            The initial conditions arrays.
        params: MappingParameters
            The integration parameters.
        """

    def run(
        self,
        qfactor: Qfactor,
        currents: Currents,
        bfield: Bfield,
        perturbation: Perturbation,
    ):
        """Integrates the particle, storing its evolution.

        Parameters
        ----------
        qfactor: Qfactor
            The equilibrium's qfactor.
        currents: Currents
            The equilibrium's plasma current.
        bfield: Bfield
            The equilibrium's magnetic field.
        per: Perturbation
            The equilibrium's perturbation.
        """
