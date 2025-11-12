#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "xarray",
#   "netcdf4",
#   "numpy",
# ]
# ///

"""Converts a `npz` file to a netcdf4 file, with the expected fields and variable names."""

import sys
from numpy._typing import _128Bit
import xarray as xr
import numpy as np
from pathlib import Path
from warnings import warn

xr.set_options(
    display_max_rows=30,
    display_width=100,
    netcdf_engine_order=["netcdf4"],
)

INPUT = Path(sys.argv[1])
OUTPUT = Path(sys.argv[2])
npz = np.load(INPUT)

print(f"Loaded .npz file from '{INPUT.absolute()}'")

# NOTE: Some arrays seem to have endianness explicitly stored in the dtype
# ('<' or '>'). Convert to the native endianness ('=', default) and expected
# type, just in case.

# NOTE: We need to divide both fluxes by 2π.

baxis = npz.get("BB")[0, 0].astype("f8")  # [Tesla]
raxis = npz.get("Rmaj").astype("f8")  # [meters]
psip = npz.get("psipol").astype("f8") / (2 * np.pi)
psi = npz.get("psitor").astype("f8") / (2 * np.pi)
theta = npz.get("theta").astype("f8")
i = npz.get("I").astype("f8")
g = npz.get("g").astype("f8")
q = npz.get("q").astype("f8")
b = npz.get("BB").astype("f8")
r = npz.get("RR").astype("f8")
z = npz.get("ZZ").astype("f8")

# If no perturbations are found in the npz file, this will create empty
# Variables, which we can drop at the end.
default_coord = np.array([])
default_array = np.full((0, 0, len(psip)), np.nan)
m = npz.get("m", default_coord).astype("i8")
n = npz.get("n", default_coord).astype("i8")
alphas = npz.get("alphas", default_array).astype("f8")
phases = npz.get("phases", default_array).astype("f8")

print("Exctracted all variables from npz file")

# Normalize
mp = 1.672621923e-27
qp = 1.602176634e-19

w0 = qp / mp * baxis  # s^-1

psip_norm = psip * mp * w0 * raxis**2 / qp
psi_norm = psi * mp * w0 * raxis**2 / qp
g_norm = g / (baxis * raxis)
i_norm = i / (baxis * raxis)
b_norm = b / baxis
alphas_norm = alphas / raxis


# ======================= Scalars


baxis_var = xr.Variable(
    data=baxis,
    dims=[],
    attrs=dict(
        description="The magnetic field strength on the axis `B0`",
        units="[T]",
    ),
)

raxis_var = xr.Variable(
    data=raxis,
    dims=[],
    attrs=dict(
        description="The major radius `R`",
        units="[m]",
    ),
)


# ======================= Coordinates

# Apart from `θ`, `ψp` and `ψ`, we also use the mode numbers `m` and `n` as
# coordinates. This way, the `a{m,n}(ψp)` 1D arrays can be easily accessed
# with `dataset.alphas[m, n]`

theta_coord = xr.Variable(
    data=theta,
    dims=["theta"],
    attrs=dict(description="Boozer theta coordinate", units="[rads]"),
)

psip_coord = xr.Variable(
    data=psip_norm,
    dims=["psip"],
    attrs=dict(description="Poloidal flux coordinate", units=["Normalized"]),
)

psi_coord = xr.Variable(
    data=psi_norm,
    dims=["psi"],
    attrs=dict(description="Toroidal flux coordinate", units=["Normalized"]),
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


# ======================= SI Variables


q_var = xr.Variable(
    data=q,
    dims=["psip"],
    attrs=dict(description="Magnetic q-factor", units="dimensionless"),
)

g_var = xr.Variable(
    data=g,
    dims=["psip"],
    attrs=dict(description="Toroidal current", units="[Tm]"),
)

i_var = xr.Variable(
    data=i,
    dims=["psip"],
    attrs=dict(description="Poloidal current", units="[Tm]"),
)

b_var = xr.Variable(
    data=b,
    dims=["psip", "theta"],
    attrs=dict(description="Magnetic field strength", units="[T]"),
)

r_var = xr.Variable(
    data=r,
    dims=["psip", "theta"],
    attrs=dict(description="Lab horizontal coordinate", units="[m]"),
)

z_var = xr.Variable(
    data=z,
    dims=["psip", "theta"],
    attrs=dict(description="Lab vertical coordinate", units="[m]"),
)


# ======================= Normalized Variables


g_norm_var = xr.Variable(
    data=g_norm,
    dims=["psip"],
    attrs=dict(description="Toroidal current normalized", units="Normalized"),
)

i_norm_var = xr.Variable(
    data=i_norm,
    dims=["psip"],
    attrs=dict(description="Poloidal current normalized", units="Normalized"),
)

b_norm_var = xr.Variable(
    data=b_norm,
    dims=["psip", "theta"],
    attrs=dict(description="Magnetic field strength", units="Normalized"),
)


# ======================= Perturbations


alphas_var = xr.Variable(
    data=alphas,
    dims=["m", "n", "psip"],
    attrs=dict(description="Mode amplitudes α{m,n}(ψp)", units="[m]"),
)

phases_var = xr.Variable(
    data=phases,
    dims=["m", "n", "psip"],
    attrs=dict(description="Mode phases φ{m,n}(ψp)", units="[rads]"),
)

alphas_norm_var = xr.Variable(
    data=alphas_norm,
    dims=["m", "n", "psip"],
    attrs=dict(description="Mode amplitudes α{m,n}(ψp)", units="Normalized"),
)


# WARN: Replace `inf` and `nan` values that may appear with 0.
np.nan_to_num(alphas_var.to_numpy(), nan=0, posinf=0, neginf=0, copy=False)
np.nan_to_num(alphas_norm_var.to_numpy(), nan=0, posinf=0, neginf=0, copy=False)


# ========================================================


COORDS = {
    "psip": psip_coord,
    "theta": theta_coord,
    "m": m_coord,
    "n": n_coord,
    "psi": psi_coord,
}

VARIABLES = {
    "baxis": baxis_var,
    "raxis": raxis_var,
    "q": q_var,
    "g": g_var,
    "i": i_var,
    "b": b_var,
    "R": r_var,
    "Z": z_var,
    "alphas": alphas_var,
    "phases": phases_var,
    "g_norm": g_norm_var,
    "i_norm": i_norm_var,
    "b_norm": b_norm_var,
    "alphas_norm": alphas_norm_var,
}

dataset = xr.Dataset(
    data_vars=VARIABLES,
    coords=COORDS,
)

# Remove empty Variables in case of no perturbations
dataset = dataset.drop_vars(
    [
        dataset[item].name
        for item in dataset.coords | dataset.data_vars
        if dataset[item].size == 0
    ]
)

print("Created dataset")

# makes `ncdump -h` cleaner
for item in dataset.coords | dataset.data_vars:
    dataset[item].attrs["_FillValue"] = False


for item in dataset.coords | dataset.data_vars:
    if not np.all(np.isfinite(dataset[item].to_numpy())):
        warn(f"NaN or inf found in `{dataset[item].name}`")

dataset.to_netcdf(OUTPUT)

print(f"Stored dataset at '{OUTPUT.absolute()}'")
