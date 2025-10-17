"""Converts a `npz` file to a netcdf4 file, with the expected fields and variable names"""

import xarray as xr
import numpy as np

# Extract arrays
npz = np.load("./usr/30840_103_axi_rev.npz")

psi_p_dataSI = npz["psipol"] / (2 * np.pi)
psi_dataSI = npz["psitor"] / (2 * np.pi)
boozer_theta_data = npz["theta"]
baxis = npz["BB"][0, 0]
raxis = npz["Rmaj"]
i = npz["I"]
g = npz["g"]
q = npz["q"]
b = npz["BB"]
s = npz["s"]
r = npz["RR"]
z = npz["ZZ"]

# Normalize
mp = 1.672621923e-27
qp = 1.602176634e-19

w0 = qp / mp * baxis  # s^-1

psi_p_data = psi_p_dataSI * mp * w0 * raxis**2 / qp
psi_data = psi_dataSI * mp * w0 * raxis**2 / qp
g_norm = g / (baxis * raxis)
i_norm = i / (baxis * raxis)
b_norm = b / baxis


# ======================= Coordinates ========================


boozer_theta_coord = xr.DataArray(
    data=boozer_theta_data,
    dims=["boozer_theta"],
    coords=dict(boozer_theta=(["boozer_theta"], boozer_theta_data)),
    attrs=dict(description="Boozer theta coordinate", units="[rads]"),
)

psi_p_coord = xr.DataArray(
    data=psi_p_data,
    dims=["psi_p"],
    coords=dict(psi_p=(["psi_p"], psi_p_data)),
    attrs=dict(description="Poloidal flux coordinate", units=["Normalized"]),
)

psi_coord = xr.DataArray(
    data=psi_data,
    dims=["psi"],
    coords=dict(psi=(["psi"], psi_data)),
    attrs=dict(description="Toroidal flux coordinate", units=["Normalized"]),
)

# ===================== 1D SI variables ======================


g_da = xr.DataArray(
    data=g,
    dims=["psi_p"],
    coords=dict(psi_p=(["psi_p"], psi_p_data)),
    attrs=dict(description="Toroidal current", units="[T][m]"),
)

i_da = xr.DataArray(
    data=i,
    dims=["psi_p"],
    coords=dict(psi_p=(["psi_p"], psi_p_data)),
    attrs=dict(description="Poloidal current", units="[T][m]"),
)

q_da = xr.DataArray(
    data=q,
    dims=["psi_p"],
    coords=dict(psi_p=(["psi_p"], psi_p_data)),
    attrs=dict(description="qfactor", units="dimensionless"),
)
# ===================== 2D SI variables ======================

b_da = xr.DataArray(
    data=b,
    dims=["psi_p", "boozer_theta"],
    coords=dict(
        psi_p=(["psi_p"], psi_p_data),
        boozer_theta=(["boozer_theta"], boozer_theta_data),
    ),
    attrs=dict(description="Magnetic field strength", units="[T]"),
)

r_da = xr.DataArray(
    data=r,
    dims=["psi_p", "boozer_theta"],
    coords=dict(
        psi_p=(["psi_p"], psi_p_data),
        boozer_theta=(["boozer_theta"], boozer_theta_data),
    ),
    attrs=dict(description="Lab horizontal coordinate", units="[m]"),
)

z_da = xr.DataArray(
    data=z,
    dims=["psi_p", "boozer_theta"],
    coords=dict(
        psi_p=(["psi_p"], psi_p_data),
        boozer_theta=(["boozer_theta"], boozer_theta_data),
    ),
    attrs=dict(description="Lab vertical coordinate", units="[m]"),
)


# ================== 1D Normalized variables ===================


g_norm_da = xr.DataArray(
    data=g_norm,
    dims=["psi_p"],
    coords=dict(psi_p=(["psi_p"], psi_p_data)),
    attrs=dict(description="Toroidal current normalized", units="Normalized"),
)

i_norm_da = xr.DataArray(
    data=i_norm,
    dims=["psi_p"],
    coords=dict(psi_p=(["psi_p"], psi_p_data)),
    attrs=dict(description="Poloidal current normalized", units="Normalized"),
)

b_field_norm_da = xr.DataArray(
    data=b_norm,
    dims=["psi_p", "boozer_theta"],
    coords=dict(
        psi_p=(["psi_p"], psi_p_data),
        boozer_theta=(["boozer_theta"], boozer_theta_data),
    ),
    attrs=dict(description="Magnetic field strength", units="Normalized"),
)


# ================== 1D Normalized variables ===================


dataset = xr.Dataset(
    data_vars=dict(
        g=g_da,
        i=i_da,
        q=q_da,
        b=b_da,
        g_norm=g_norm_da,
        I_norm=i_norm_da,
        b_field_norm=b_field_norm_da,
        R=r_da,
        Z=z_da,
        Baxis=baxis,
        raxis=raxis,
        zaxis=np.nan,
    ),
    coords=dict(
        psi_p=psi_p_coord,
        boozer_theta=boozer_theta_coord,
        psi=psi_coord,
    ),
)

for var in dataset.variables:
    dataset[var].attrs["_FillValue"] = False

dataset.to_netcdf(
    "data.nc",
    engine="netcdf4",
)
