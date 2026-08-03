# NetCDF Convention

Version 0.0.1

DEXTER reads equilibrium data from a [netCDF] file. The variables must follow the following conventions:

## Scalars

* `baxis`: The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
* `raxis`: The horizontal position of the magnetic axis $R_0$ $[m]$.
* `zaxis`: The horizontal position of the magnetic axis $[m]$.
* `rgeo`: The geometrical axis (device major radius) in $[m]$.

!!! note "Normalizations"

    For the normalizations, `baxis` and `raxis` should be used. They correspond to the center of the smallest flux surface. `rgeo` is characteristic of the device, while `raxis` depends on the configuration.

## Coordinates

!!! tip

    If either `psi_norm` or `psip_norm` is missing from the netCDF file, it is calculated from the other by integrating $q(\psi_p)$ or $\iota(\psi)$ respectively. In the case that the calculated values are monotonic, the other flux can be used as a flux coordinate as well.

* `psi_norm`: The toroidal flux coordinate $\psi$ divided by $2\pi$, in Normalized Units (consistent with the `eqdisk` files).
* `psip_norm`: The poloidal flux coordinate $\psi_p$ divided by $2\pi$, in Normalized Units (consistent with the `eqdisk` files).
* `theta`: The Boozer theta coordinate $\theta$ in $[rads]$.
* `m`: The poloidal mode number $m$ (*index* coordinate).
* `n`: The toroidal mode number $n$ (*index* coordinate).

Not used in any calculations:

* `r_norm`: The radial distance coordinate $r$ in Normalized Units.

## Variables

!!! note "Magnetic flux coordinate"

    Here, $\psi$ and $\psi_p$ are used interchangeably.

* `q`: The safety factor $q(\psi)$.
* `g_norm`: The toroidal plasma current $g(\psi)$ in Normalized Units.
* `i_norm`: The poloidal plasma current $I(\psi)$ in Normalized Units.
* `b_norm`: The magnetic field strength $B(\psi, \theta)$ in Normalized Units.
* `alphas_norm`: The harmonic amplitudes $\alpha_{m,n}(\psi)$ in Normalized Units.
* `phases`: The harmonic phases $\phi_{m,n}(\psi)$ in $[rads]$.

Original SI data (not used in any calculations):

* `r`: The radial distance coordinate $r(\psi)$ in $[m]$.
* `g`: The toroidal plasma current $g(\psi)$ in $[T \cdot m]$.
* `I`: The poloidal plasma current $I(\psi)$ in $[T \cdot m]$.
* `b`: The magnetic field strength $B(\psi, \theta)$ in $[T]$.
* `jacobian`: The Jacobian matrix $J(\psi, \theta)$ in $[m/T]$.
* `alphas`: The harmonic amplitudes $\alpha_{m,n}(\psi)$ in $[m]$.
* `rlab`: The lab horizontal coordinate $R(\psi, \theta)$ in $[m]$.
* `zlab`: The lab vertical coordinate $Z(\psi, \theta)$ in $[m]$.
* `psip`: The poloidal flux coordinate $\psi_p$ in $[T \cdot m^2]$.
* `psi`: The toroidal flux coordinate $\psi$ in $[T \cdot m^2]$.

## Attributes

Optional file attributes to avoid confusion.

* `description`: A short description.
* `date`: The file's creation date.
* `script`: The script used to create the file.
* `version`: The *convention* version, using [Semantic Versioning].


[netCDF]: https://www.unidata.ucar.edu/software/netcdf
[Semantic Versioning]: https://semver.org/
