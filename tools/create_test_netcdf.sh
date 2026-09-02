# Creates all testing `netCDF` files and places them at the correct location.

python ./tools/create_test_netcdf.py -o ./crates/dexter-equilibrium/test_netcdf.nc -c both
python ./tools/create_test_netcdf.py -o ./crates/dexter-equilibrium/toroidal_test_netcdf.nc -c toroidal
python ./tools/create_test_netcdf.py -o ./crates/dexter-equilibrium/poloidal_test_netcdf.nc -c poloidal
ln -srvf ./crates/dexter-equilibrium/test_netcdf.nc ./crates/dexter-equilibrium/netcdf.nc

ln -srvf ./crates/dexter-equilibrium/test_netcdf.nc ./netcdf.nc
ln -srvf ./crates/dexter-equilibrium/test_netcdf.nc ./crates/dexter-simulate/test_netcdf.nc
ln -srvf ./crates/dexter-equilibrium/test_netcdf.nc ./crates/dexter-simulate/netcdf.nc
ln -srvf ./crates/dexter-equilibrium/toroidal_test_netcdf.nc ./crates/dexter-simulate/toroidal_test_netcdf.nc
ln -srvf ./crates/dexter-equilibrium/poloidal_test_netcdf.nc ./crates/dexter-simulate/poloidal_test_netcdf.nc
