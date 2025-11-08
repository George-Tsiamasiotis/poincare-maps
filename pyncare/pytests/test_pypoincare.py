from pyncare import Qfactor, Currents, Bfield, Perturbation, Poincare
import numpy as np


def test_poincare_fields(poincare: Poincare):
    assert isinstance(poincare.alpha, float)
    assert isinstance(poincare.section, str)
    assert poincare.section in ["ConstTheta", "ConstZeta"]
    assert isinstance(poincare.angles, np.ndarray)
    assert isinstance(poincare.fluxes, np.ndarray)


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
