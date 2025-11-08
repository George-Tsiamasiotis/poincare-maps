from pyncare import Harmonic, Perturbation


def test_pypertyrbation_eval(perturbation: Perturbation):
    psip = 0.015
    theta = 1
    zeta = 2
    assert isinstance(perturbation.p(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dpsip(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dtheta(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dzeta(psip, theta, zeta), float)
    assert isinstance(perturbation.dp_dt(psip, theta, zeta), float)


def test_perturbation_getitem(perturbation: Perturbation):
    assert isinstance(perturbation[0], Harmonic)
    assert isinstance(perturbation[1], Harmonic)


def test_repr(perturbation: Perturbation):
    str(perturbation)
