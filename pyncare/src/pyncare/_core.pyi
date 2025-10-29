import numpy as np

class Qfactor:
    """q-factor reconstructed from a NetCDF file.

    Attributes
    ----------
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psi_wall: float
        The value of the toroidal flux ψ at the wall.
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    """

    def __init__(self, path: str, typ: str):
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

    def psip_data(self) -> np.ndarray:
        """The NetCDF ψp data used to construct the q(ψp) and ψ(ψp) splines."""

    def q_data(self) -> np.ndarray:
        """The NetCDF q data used to construct the q(ψp) spline."""

    def psi_data(self) -> np.ndarray:
        """The NetCDF ψp data used to construct the ψ(ψp) spline."""

    def q_data_derived(self) -> np.ndarray:
        """The q values, as calculated from dψ/dψp, at the ψp data."""

class Current:
    """Plasma current reconstructed from a NetCDF file.

    Attributes
    ----------
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    """

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

    def psip_data(self) -> np.ndarray:
        """The NetCDF ψp data used to construct the g(ψp) and I(ψp) splines."""

    def g_data(self) -> np.ndarray:
        """The NetCDF g data used to construct the g(ψp) spline."""

    def i_data(self) -> np.ndarray:
        """The NetCDF I data used to construct the I(ψp) spline."""

    def dg_dpsip_data(self) -> np.ndarray:
        """The dg/dψp values evaluated at the ψp data"""

    def di_dpsip_data(self) -> np.ndarray:
        """The dI/dψp values evaluated at the ψp data"""

class Bfield:
    """Magnetic field reconstructed from a NetCDF file.

    Attributes
    ----------
    baxis: float
        The magnetic field strength on the axis in [T].
    raxis: float
        The major radius in [m].
    psip_wall: float
        The value of the poloidal flux ψp at the wall.
    psi_wall: float
        The value of the toroidal flux ψ at the wall.
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation.
    """

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

    def psip_data(self) -> np.ndarray:
        """The NetCDF ψp data used to construct the b(ψp, θ) spline."""

    def theta_data(self) -> np.ndarray:
        """The NetCDF θ data used to construct the b(ψp, θ) spline."""

    def b_data(self) -> np.ndarray:
        """The NetCDF b data used to construct the b(ψp, θ) spline."""

    def r_data(self) -> np.ndarray:
        """The NetCDF R data used to construct the R(ψp, θ) spline."""

    def z_data(self) -> np.ndarray:
        """The NetCDF Z data used to construct the Z(ψp, θ) spline."""

    def db_dpsip_data(self) -> np.ndarray:
        """The db/dψp values evaluated at the (ψp, θ) data"""

    def db_dtheta_data(self) -> np.ndarray:
        """The db/dψp values evaluated at the (ψp, θ) data"""
