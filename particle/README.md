## Particle
Defines the `Particle` object, a representation of a particle in a toroidal fusion device, in Normalized Units space:

- mass = charge = 1
- All quantities are in Normalized Units

The Particle can be integrated with Runge-Kutta-Fehlberg 4(5) routine, with 2 methods of calculating the error and next optimal step size (see [features](#Features))

### Integration

The time evolution of the particle is calculated by straight forward integration from a set of initial conditions.

### Mapping

The particle is integrated in the same way, but only its intersections with the constant surfaces α=θ or α=ζ are stored.
The intersections are calculated using [`Hénon's trick`] which guarantees an error smaller than the solver's tolerance.

### Features

- `phase-interpolation`: See [`phase-interpolation`].
- `energy-adaptive-step`: The solver's step-size is adjusted by trying to keep the Energy difference between every step under a certain threshold. Otherwise, it is adjusted by estimating the local truncation error (default, faster).Both tolerances can be changed in the configuration crate.

[`Hénon's trick`]:https://doi.org/10.1016/0167-2789(82)90034-3
[`phase-interpolation`]: ../equilibrium/index.html#features