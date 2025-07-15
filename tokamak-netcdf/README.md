## tokamak-netcdf

A package for opening [`NetCDF`] files of reconstructed [`Tokamak`] equilibria and providing a more suitable interface for accessing the data and performing calculations.


## Notes

Requires the [`netCDF-C`] library, which is available in most linux package managers. In case it is not, it can be statically linked with the 'static' feature, which is provided by the [`netcdf crate`].

[`NetCDF`]: https://www.unidata.ucar.edu/software/netcdf
[`netCDF-C`]: https://github.com/Unidata/netcdf-c
[`netcdf crate`]: https://github.com/georust/netcdf
[`Tokamak`]: https://en.wikipedia.org/wiki/Tokamak