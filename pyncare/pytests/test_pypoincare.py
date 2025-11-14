from pyncare import Qfactor, Currents, Bfield, Perturbation, Poincare
import numpy as np

from pyncare._core import Particle


def test_poincare_fields(poincare: Poincare):
    assert isinstance(poincare.alpha, float)
    assert isinstance(poincare.section, str)
    assert poincare.section in ["ConstTheta", "ConstZeta"]
    assert isinstance(poincare.angles, np.ndarray)
    assert isinstance(poincare.fluxes, np.ndarray)
    assert isinstance(poincare.zetas, np.ndarray)
    assert isinstance(poincare.psips, np.ndarray)
    assert isinstance(poincare.thetas, np.ndarray)
    assert isinstance(poincare.psis, np.ndarray)
    assert isinstance(poincare[1], Particle)


def test_pypoincare_run(
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
    perturbation: Perturbation,
    poincare: Poincare,
):
    poincare.run(
        qfactor=qfactor,
        currents=currents,
        bfield=bfield,
        perturbation=perturbation,
    )


def test_repr(poincare: Poincare):
    str(poincare)
