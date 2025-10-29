import numpy as np

class Qfactor:
    """
    q-factor reconstructed from a NetCDF file.

    Parameters
    ----------
    path: str
        The path to the NetCDF file.
    typ: str
        The type of Interpolation. Available types: "Linear", "Cubic", "Akima", "AkimaPeriodic",
        "Steffen".

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
        """Constructor"""

    def q(psip: float) -> float:
        """The q value evaluated at ψp"""

    def psi(psip: float) -> float:
        """The ψ value evaluated at ψp"""

    def psip_data() -> np.ndarray:
        """The NetCDF ψp data used to construct the q(ψp) and ψ(ψp) splines."""

    def q_data() -> np.ndarray:
        """The NetCDF q data used to construct the q(ψp) spline."""

    def psi_data() -> np.ndarray:
        """The NetCDF ψp data used to construct the ψ(ψp) spline."""

    def q_data_derived() -> np.ndarray:
        """The q values, as calculated from dψ/dψp, at the ψp data."""
