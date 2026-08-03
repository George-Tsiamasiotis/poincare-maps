# DEXTER

___Dynamics of Experimental Toroidal Equilibrium Reconstructions___


A code for performing calculations and simulations of particles and magnetic field lines, with reconstructed equilibria from experimental data.

The data must be in [netCDF] format.

The bulk computations are implemented in *[Rust]*. The Rust code consists of standalone *[crates]* that can be used independently and as is.

The Python interface directly exposes all underlying objects and routines in the form of a single python package:

```python
>>> # TODO:
```

while also providing plotting methods and scripts for handling and converting netCDF files.

[netCDF]: https://www.unidata.ucar.edu/software/netcdf
[Rust]: https://rust-lang.org/
[crates]: https://rustup.dev/fundamentals/crates.html#crates
