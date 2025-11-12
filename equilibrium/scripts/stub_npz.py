#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "numpy",
# ]
# ///

"""Creates a stub .npz file, in the expected format.

Output file can be converted to a stub NetCDF file with `./equilibrium/scripts/npz_to_netcdf.py`.
"""

import sys
import numpy as np
from pathlib import Path

OUTPUT = Path(sys.argv[1])

baxis = 1.5
Rmaj = 2
theta = np.linspace(0, np.pi, 200)
psipol = np.linspace(0, 3, 100)
psitor = np.linspace(0, 1, 100)
I = np.linspace(0, 2, 100)
g = np.linspace(2, 0, 100)
q = np.linspace(1, 2, 100)
BB = np.random.random((len(psipol), len(theta)))
RR = np.random.random((len(psipol), len(theta)))
ZZ = np.random.random((len(psipol), len(theta)))
empty = np.array([])

m = np.arange(-1, 4, dtype="i8")
n = np.arange(-2, 8, dtype="i8")
alphas = np.random.random((len(m), len(n), len(psipol)))
phases = np.random.random((len(m), len(n), len(psipol)))
# Inject the (2, 3) mode with specific values to test that the extraction is done correctly.
# If either the indexes or the values change, the test must be updated accordingly.
alphas[2, 3, 0] = 1111
phases[2, 3, 0] = 9999
alphas[2, 3, -1] = 11111
phases[2, 3, -1] = 99999


np.savez(
    file=OUTPUT,
    baxis=baxis,
    Rmaj=Rmaj,
    theta=theta,
    psipol=psipol,
    psitor=psitor,
    I=I,
    g=g,
    q=q,
    BB=BB,
    RR=RR,
    ZZ=ZZ,
    m=m,
    n=n,
    alphas=alphas,
    phases=phases,
)

print(f"Created stub NetCDF file at '{OUTPUT.absolute()}'")
