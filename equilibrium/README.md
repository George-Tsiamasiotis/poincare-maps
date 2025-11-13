## Equilibrium

Handles the data extraction from [`NetCDF`] files, as well as the reconstruction of the equilibrium, providing interpolation methods for calculating all the relevant quantities.

This crate requires the [`netCDF-C`] library, which is available in most linux package managers.

`libnetcdf` can be statically linked with the `netcdf-static` feature, which is provided by the
[`netcdf crate`].

### NetCDF files

`equilibrium/scripts/npz_to_netcdf.py` converts a `.npz` file to a [`NetCDF`] file, adding the normalized data, adding `description` and `units` attributes in every variable and filling possible missing values. It should be used to create a stub .nc file for testing:

### Features

- `phase-interpolation`: Calculates the phase α(ψp) at each step, by interpolating over the 1D phase arrays in the NetCDF file. Otherwise, the average value of the array is calculated beforehand and used as a constant phase throughout the integration.  

### Testing

`equilibrium/scripts/stub_npz.py` creates a stub `.npz` file, containing the expected arrays in their expected names and shapes.

The stub .nc file can then be used for testing the crate, as it has the same layout as a real one:

```bash
// from project root
./equilibrium/scripts/stub_npz.py "./data/stub_npz.npz" 
./equilibrium/scripts/npz_to_netcdf.py "./data/stub_npz.npz" "./data/stub_netcdf.nc"
cargo test -p equilibrium
```


[`netCDF`]: https://www.unidata.ucar.edu/software/netcdf
[`netCDF-C`]: https://github.com/Unidata/netcdf-c
[`netcdf crate`]: https://github.com/georust/netcdf
