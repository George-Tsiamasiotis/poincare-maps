# Geometry Objects

::: dexter.GeometryObject

---

::: dexter.LarGeometry
    options:
        show_bases: true
        inherited_members: false

A LAR configuration's geometry is described by the following formulas:

$$
r(\psi) = \sqrt{2\psi}
$$

$$
\psi(r) = r^2/2
$$

$$
R(\psi, \theta) = R_{geo} + R_{geo}\sqrt{2\psi}\cos\theta
$$

$$
Z(\psi, \theta) = R_{geo}\sqrt{2\psi}\sin\theta
$$

---

::: dexter.NcGeometry
    options:
        show_bases: true
        inherited_members: false
