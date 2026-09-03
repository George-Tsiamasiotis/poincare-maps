# Creates all testing `netCDF` files and places them at the correct location.

python ./tools/create_test_netcdf.py -o ./crates/dexter-machine/test_netcdf.nc -c both
python ./tools/create_test_netcdf.py -o ./crates/dexter-machine/toroidal_test_netcdf.nc -c toroidal
python ./tools/create_test_netcdf.py -o ./crates/dexter-machine/poloidal_test_netcdf.nc -c poloidal
ln -srvf ./crates/dexter-machine/test_netcdf.nc ./crates/dexter-machine/netcdf.nc

ln -srvf ./crates/dexter-machine/test_netcdf.nc ./netcdf.nc
ln -srvf ./crates/dexter-machine/test_netcdf.nc ./crates/dexter-simulate/test_netcdf.nc
ln -srvf ./crates/dexter-machine/test_netcdf.nc ./crates/dexter-simulate/netcdf.nc
ln -srvf ./crates/dexter-machine/toroidal_test_netcdf.nc ./crates/dexter-simulate/toroidal_test_netcdf.nc
ln -srvf ./crates/dexter-machine/poloidal_test_netcdf.nc ./crates/dexter-simulate/poloidal_test_netcdf.nc
