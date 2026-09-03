# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "semver",
#     "numpy",
#     "xarray",
# ]
# ///

"""Creates a test netCDF file with a LAR configuration.

Flux coordinates behavior:

    1. -c = "both": equispaced ψ values, parabolic q.
    2. -c = "toroidal": equispaced ψ values, ψp = sin(2πψ), and q(ψ) = dψp/dψ.
    3. -c = "poloidal": equispaced ψp values, ψ = sin(2πψp), and q(ψ) = dψp/dψ = 1/(dψ/dψ).

"""

import argparse
import numpy as np
import xarray as xr
import tomllib
from math import sqrt
from semver import Version
from pathlib import Path
from datetime import datetime
from typing import assert_never

with open("./Cargo.toml", "rb") as f:
    VERSION = Version.parse(tomllib.load(f)["workspace"]["metadata"]["netcdf_version"])

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"netCDF version {VERSION}",
)
parser.add_argument(
    "-o",
    "--output",
    help="the output file name (leave empty when running interactively)",
    type=str,
    default="",  # when run interactively
)
parser.add_argument(
    "-c",
    "--coord",
    help="the good (monotonic) flux coordinate(s)",
    type=str,
    default="both",
    choices=["both", "toroidal", "poloidal"],
)
args = parser.parse_args()

# ==========================================================================

OUTPUT = Path(args.output)

#########
# Scalars
#########

baxis = 1.5
raxis = 1.75
zaxis = 0
rgeo = raxis
a = 0.5  # "small radius"

flux_last_value_norm = 0.5  # can be either ψ_last or ψp_last

##################################
# Coordinates, Fluxes and q-factor
##################################

FLUX_SURFACES = 100
THETAS = 200

theta = np.linspace(0, 2 * np.pi, THETAS)
m = np.asarray([2, 3])
n = np.asarray([1, 2])

match args.coord:
    case "both":  # Parabolic
        qaxis = 1.1
        qlast = 3.9
        psi_norm = np.linspace(0, flux_last_value_norm, FLUX_SURFACES)
        atan_arg = psi_norm * sqrt(qlast - qaxis) / (flux_last_value_norm * sqrt(qaxis))
        coef = flux_last_value_norm / sqrt(qaxis * (qlast - qaxis))
        psip_norm = coef * np.atan(atan_arg)
        q = qaxis + (qlast - qaxis) * (psi_norm / flux_last_value_norm) ** 2
    case "toroidal":  # equispaced ψ values, ψp = sin(2πψ), and q(ψ) = dψp/dψ.
        psi_norm = np.linspace(0, flux_last_value_norm, FLUX_SURFACES)
        psip_norm = np.sin(2 * np.pi * psi_norm)  # no longer monotonic in [0, psi_last]
        q = np.gradient(psip_norm, psi_norm)
        assert not np.all(np.diff(psip_norm) > 0)
    case "poloidal":  # equispaced ψp values, ψ = sin(2πψp), and q(ψ) = 1/(dψp/dψ).
        psip_norm = np.linspace(0, flux_last_value_norm, FLUX_SURFACES)
        psi_norm = 1.0 - np.sin(
            2 * np.pi * psip_norm
        )  # no longer monotonic in [0, psi_last]
        q = np.gradient(psi_norm, psip_norm)  # inverse
        assert not np.all(np.diff(psi_norm) > 0)
    case _:
        assert_never(args.coord)

# `r` is somewhat artificial, since we cannot always calculate it from the fluxes
r_norm = np.linspace(0, a, FLUX_SURFACES) / raxis


# Currents
g_norm = np.ones(FLUX_SURFACES)
i_norm = np.zeros(FLUX_SURFACES)

# LAR Bfield
r_norm_grid, theta_grid = np.meshgrid(r_norm, theta, indexing="ij")
b_norm = 1 - r_norm_grid * np.cos(theta_grid)

# Lab coordinates
rlab = raxis + r_norm_grid * np.cos(theta_grid)
zlab = r_norm_grid * np.sin(theta_grid)
jacobian = ((g_norm.T + i_norm.T / q.T) / b_norm.T**2).T  # Unsure

###############
# Perturbations
###############

# Should be just small variations
a12 = a13 = a22 = a23 = 1e-3 * np.cos(2 * psip_norm)
p12 = p13 = p22 = p23 = np.sin(10 * psip_norm)
alphas_norm = np.asarray([[a12, a23], [a13, a23]])
phases = np.asarray([[p12, p23], [p13, p23]])

# Value Injection
# Inject the (2, 2) mode with specific values to test that the extraction
# is done correctly. If either the indexes or the values change, the test
# must be updated accordingly.
phases[0, 1, 0] = 9999
phases[0, 1, -1] = 99999
alphas_norm[0, 1, 0] = 1111
alphas_norm[0, 1, -1] = 11111


#################
# Misc quantities
#################

psi = psi_norm * (baxis * raxis**2)
psip = psip_norm * (baxis * raxis**2)
r = r_norm * raxis
g = g_norm * baxis * raxis
i = i_norm * baxis * raxis
b = b_norm / baxis
alphas = alphas_norm * raxis

# ======================= Scalars

baxis_var = xr.Variable(
    data=baxis,
    dims=[],
    attrs=dict(
        description="The magnetic field strength on the axis `B0`",
        units="[T]",
        normalization="This value should be used for normalizations",
    ),
)

raxis_var = xr.Variable(
    data=raxis,
    dims=[],
    attrs=dict(
        description="The horizontal position of the magnetic axis R0",
        units="[m]",
        normalization="This value should be used for normalizations",
    ),
)

zaxis_var = xr.Variable(
    data=zaxis,
    dims=[],
    attrs=dict(description="The vertical position of the magnetic axis", units="[m]"),
)

rgeo_var = xr.Variable(
    data=rgeo,
    dims=[],
    attrs=dict(description="The geometrical axis (device major radius)", units="[m]"),
)

# ======================= Coordinates

theta_coord = xr.Variable(
    data=theta,
    dims=["theta"],
    attrs=dict(description="Boozer theta coordinate", units="[rads]"),
)

psi_norm_coord = xr.Variable(
    data=psi_norm,
    dims=["psi_norm"],
    attrs=dict(description="Toroidal flux coordinate", units="Normalized"),
)

m_coord = xr.Variable(
    data=m,
    dims=["m"],
    attrs=dict(description="Poloidal mode number"),
)

n_coord = xr.Variable(
    data=n,
    dims=["n"],
    attrs=dict(description="Toroidal mode number"),
)

psip_norm_coord = xr.Variable(
    data=psip_norm,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal flux coordinate", units="Normalized"),
)

r_norm_coord = xr.Variable(
    data=r_norm,
    dims=["r_norm"],
    attrs=dict(description="Radial coordinate", units="Normalized"),
)


# ======================= SI Variables

psip_var = xr.Variable(
    data=psip,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal flux coordinate", units="[Tm]"),
)

r_var = xr.Variable(
    data=r,
    dims=["r_norm"],
    attrs=dict(description="Radial coordinate", units="[m]"),
)

q_var = xr.Variable(
    data=q,
    dims=["psip_norm"],
    attrs=dict(description="Magnetic q-factor", units="dimensionless"),
)

psi_var = xr.Variable(
    data=psi,
    dims=["psi_norm"],
    attrs=dict(description="Toroidal flux coordinate", units="[Tm]"),
)


g_var = xr.Variable(
    data=g,
    dims=["psip_norm"],
    attrs=dict(description="Toroidal current", units="[Tm]"),
)

i_var = xr.Variable(
    data=i,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal current", units="[Tm]"),
)

b_var = xr.Variable(
    data=b,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Magnetic field strength", units="[T]"),
)

jacobian_var = xr.Variable(
    data=jacobian,
    dims=["psip_norm", "theta"],
    attrs=dict(description="VMEC to Boozer Jacobian", units="[m/T]"),
)

rlab_var = xr.Variable(
    data=rlab,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Lab horizontal coordinate", units="[m]"),
)

zlab_var = xr.Variable(
    data=zlab,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Lab vertical coordinate", units="[m]"),
)


# ======================= Normalized Variables


g_norm_var = xr.Variable(
    data=g_norm,
    dims=["psip_norm"],
    attrs=dict(description="Toroidal current normalized", units="Normalized"),
)

i_norm_var = xr.Variable(
    data=i_norm,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal current normalized", units="Normalized"),
)

b_norm_var = xr.Variable(
    data=b_norm,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Magnetic field strength", units="Normalized"),
)


# ======================= Perturbations


alphas_var = xr.Variable(
    data=alphas,
    dims=["m", "n", "psip_norm"],
    attrs=dict(description="Mode amplitudes α{m,n}(ψp)", units="[m]"),
)


phases_var = xr.Variable(
    data=phases,
    dims=["m", "n", "psip_norm"],
    attrs=dict(description="Mode phases φ{m,n}(ψp)", units="[rads]"),
)

alphas_norm_var = xr.Variable(
    data=alphas_norm,
    dims=["m", "n", "psip_norm"],
    attrs=dict(description="Mode amplitudes α{m,n}(ψp)", units="Normalized"),
)


# ========================================================


COORDS = {
    "psip_norm": psip_norm_coord,
    "theta": theta_coord,
    "m": m_coord,
    "n": n_coord,
    "psi_norm": psi_norm_coord,
    "r_norm": r_norm,
}

VARIABLES = {
    "baxis": baxis_var,
    "raxis": raxis_var,
    "zaxis": zaxis_var,
    "rgeo": rgeo_var,
    "psip": psip_var,
    "psi": psi_var,
    "r": r_var,
    "q": q_var,
    "g": g_var,
    "i": i_var,
    "b": b_var,
    "jacobian": jacobian_var,
    "rlab": rlab_var,
    "zlab": zlab_var,
    "alphas": alphas_var,
    "phases": phases_var,
    "g_norm": g_norm_var,
    "i_norm": i_norm_var,
    "b_norm": b_norm_var,
    "alphas_norm": alphas_norm_var,
}

match args.coord:
    case "both":
        description = "LAR testing configuration with monotonic ψ and ψp"
    case "toroidal":
        description = "LAR testing configuration with monotonic ψ and non-monotonic ψp"
    case "poloidal":
        description = "LAR testing configuration with non-monotonic ψ and monotonic ψp"
    case _:
        assert_never(args.coord)

ATTRS = {
    "description": description,
    "date": str(datetime.now().astimezone().replace(microsecond=0).isoformat()),
    "script": str(Path(__file__).stem),
    "version": str(VERSION),
}

dataset = xr.Dataset(
    data_vars=VARIABLES,
    coords=COORDS,
    attrs=ATTRS,
)


GREEN = "\033[92m"
if args.output != "":
    dataset.to_netcdf(OUTPUT)
    print(f"{GREEN}Stored dataset at '{OUTPUT.absolute()}'")
else:
    print(f"{GREEN}Created dataset (not stored)")


dataset.close()
raise SystemExit
