# Mode Objects

+ [`FluteMode`](#dexter.FluteMode): A simple analytical flute mode.
+ [`NcFluteMode`](#dexter.NcFluteMode): Single perturbation flute mode from a netCDF file.

---

::: dexter.FluteMode
    options:
        show_bases: true
        inherited_members: false

A flute mode is defined as:

$$
m(\psi, \theta, \zeta) = \epsilon\sqrt{\dfrac{\psi}{\psi_{LCFS}}}\cos(m\theta-n\zeta+\phi)
$$

---

::: dexter.NcFluteMode
    options:
        show_bases: true
        inherited_members: false

A numerical flute mode is defined as:

$$
m(\psi, \theta, \zeta) = \sum_{m,n} \alpha(\psi, \theta, \zeta)\cos\big(m\theta-n\zeta+\phi(\psi, \theta, \zeta)\big)
$$
