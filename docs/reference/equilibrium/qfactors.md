# q-factor Objects

+ [`UnityQfactor`](#dexter.UnityQfactor): Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.
+ [`ParabolicQfactor`](#dexter.ParabolicQfactor): Analytical parabolic q-factor profile.
+ [`NcQfactor`](#dexter.NcQfactor): Numerical q-factor profile reconstructed from a netCDF file.

---

::: dexter.UnityQfactor
    options:
        show_bases: true
        inherited_members: false

---

::: dexter.ParabolicQfactor
    options:
        show_bases: true
        inherited_members: false

The ParabolicQfactor's profile is described by the following formulas:

$$
q(\psi) = q_{axis} + (q_{LCFS} - q_{axis})
    \bigg( \dfrac{\psi}{\psi_{LCFS}} \bigg)^2
$$

$$
\psi_p(\psi) = \dfrac{\psi_{LCFS}}{\sqrt{q_{axis}(q_{LCFS} - q_{axis})}}
    \arctan\bigg[ \dfrac{\psi\sqrt{q_{LCFS} - q_{axis}}}{\psi_{LCFS}\sqrt{q_{axis}}} \bigg]
$$

$$
\psi(\psi_p) = \dfrac{\psi_{LCFS}\sqrt{q_{axis}}}{\sqrt{q_{LCFS} - q_{axis}}}
    \tan\bigg[ \dfrac{\sqrt{q_{axis}(q_{LCFS} - q_{axis})}}{\psi_{LCFS}}\psi_p \bigg]
$$

$$
\dfrac{d\psi_p(\psi)}{d\psi} = ... = \dfrac{1}{q(\psi)} = \iota(\psi)
$$

$$
\dfrac{d\psi(\psi_p)}{d\psi_p} =
    \dfrac{q_{axis}}{\cos^2 \bigg[
    \dfrac{\sqrt{q_{axis}(q_{LCFS} - q_{axis})}}{\psi_{LCFS}}\psi_p
    \bigg]}
    \overset{*}{=}
    q(\psi_p)
$$

$$
q(\psi_p) = q_{axis} + q_{axis} \tan^2
    \bigg[ \dfrac{\sqrt{q_{axis}(q_{LCFS}-q_{axis})}}{\psi_{LCFS}} \psi_p \bigg]
$$

$^*$ Identity: $\dfrac{1}{\cos^2\theta} = 1 + \tan^2\theta$

---

::: dexter.NcQfactor
    options:
        show_bases: true
        inherited_members: false
