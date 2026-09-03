# DEXTER - Dynamics of Experimental Toroidal Equilibrium Reconstructions

A code for performing calculations and simulations of particles and magnetic field lines, with reconstructed configurations from experimental data.

## Dependencies

1. [`libnetcdf`](https://www.unidata.ucar.edu/software/netcdf) for handling netCDF files.
2. [`openblas`](https://github.com/OpenMathLib/OpenBLAS?tab=readme-ov-file#openblas) for calculating spline coefficients.

Both are available in most Linux distributions. If compiled from source, both libraries can be linked statically with the `netcdf-static` and `openblas-static` features on the `dexter-python` crate.

## Documentation
Documentation can be found [here](https://dexter.tsiamasiotis.gr).
