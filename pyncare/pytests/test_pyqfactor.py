import pytest
from pyncare import Qfactor


def test_pyqfactor_derived_fields(qfactor: Qfactor):
    assert isinstance(qfactor.path, str)
    assert isinstance(qfactor.typ, str)
    assert isinstance(qfactor.psip_wall, float)
    assert isinstance(qfactor.psi_wall, float)


def test_pyqfactor_eval(qfactor: Qfactor):
    psip = 0.015
    assert isinstance(qfactor.q(psip), float)
    assert isinstance(qfactor.psi(psip), float)


def test_data_extraction(qfactor: Qfactor):
    assert qfactor.psip_data.ndim == 1
    assert qfactor.q_data.ndim == 1
    assert qfactor.psi_data.ndim == 1
    assert qfactor.q_data_derived.ndim == 1


def test_immutability(qfactor: Qfactor):
    with pytest.raises(AttributeError):
        qfactor.psip_wall += 1
        qfactor.psi_wall += 1
        qfactor.path = ""
        qfactor.typ = ""


def test_repr(qfactor: Qfactor):
    str(qfactor)
