# Machine

Representations of analytical and numerical machines.

Handles the data extraction from [`NetCDF`] files, as well as the reconstruction of the machine, providing interpolation methods for calculating all the relevant quantities.

This crate requires the [`netCDF-C`] library, which is available in most Linux package managers.

`libnetcdf` can be statically linked with the `netcdf-static` feature, which is provided by the
[`netcdf crate`].

## `NetCDF` files

The netCDF file must follow a [`specific convention`](https://dexter.tsiamasiotis.gr/netcdf).

[`netCDF`]: https://www.unidata.ucar.edu/software/netcdf
[`netCDF-C`]: https://github.com/Unidata/netcdf-c
[`netcdf crate`]: https://github.com/georust/netcdf
