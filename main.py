import poincare_maps as pm

bfield = pm.Bfield("./data.nc", "bicubic")
qfactor = pm.Qfactor("./data.nc", "akima")
current = pm.Current("./data.nc", "bicubic")

init = pm.InitialConditions(
    t0=0,
    theta0=0,
    psip0=0.05,
    rho0=0,
    zeta0=0,
    mu=1e-5,
    pzeta=-0.01,
)
