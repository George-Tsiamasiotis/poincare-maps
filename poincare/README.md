## Poincare

Calculates Poincare maps, by [`mapping`] particles created from many initial condtions.

[`mapping`]: ../particle/index.html#mapping

### Features

- `phase-interpolation`: See [`phase-interpolation`].
- `energy-adaptive-step`: See [`energy-adaptive-step`].


[`phase-interpolation`]: ../equilibrium/index.html#features
[`energy-adaptive-step`]: ../particle/index.html#features

### Environment variables

- `RAYON_NUM_THREADS`: The number of threads to use for the calculation. Every particle is mapped independently on its own thread.