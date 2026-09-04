import pytest
import dexter as dex
from math import isclose


def test_analytical():
    lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    machine = dex.Machine(
        geometry=dex.LarGeometry(2, 1.75, 0.5),
        qfactor=dex.UnityQfactor(lcfs),
        current=dex.LarCurrent(),
        bfield=dex.LarBfield(),
    )
    machine.perturbation = dex.Perturbation(
        [
            dex.FluteMode(1e-5, lcfs, 3, 1, 0),
            dex.FluteMode(1e-5, lcfs, 3, 2, 0),
        ]
    )
    assert machine.baxis == 2
    assert machine.raxis == 1.75
    assert machine.rlast == 0.5
    assert machine.zaxis == 0
    assert machine.rgeo == machine.raxis
    assert isinstance(machine.psi_last, float)
    assert isinstance(machine.psip_last, float)
    assert isinstance(machine.geometry, dex.GeometryObject)
    assert isinstance(machine.qfactor, dex.QfactorObject)
    assert isinstance(machine.current, dex.CurrentObject)
    assert isinstance(machine.bfield, dex.BfieldObject)
    assert isinstance(machine.perturbation, dex.Perturbation)
    machine.__repr__()
    machine.__str__()


def test_analytical_no_geometry():
    lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    machine = dex.Machine(
        qfactor=dex.UnityQfactor(lcfs),
        current=dex.LarCurrent(),
        bfield=dex.LarBfield(),
        perturbation=dex.Perturbation(
            [
                dex.FluteMode(1e-5, lcfs, 3, 1, 0),
                dex.FluteMode(1e-5, lcfs, 3, 2, 0),
            ]
        ),
    )
    assert isinstance(machine.psi_last, float)
    assert isinstance(machine.psip_last, float)
    with pytest.raises(AttributeError):
        machine.geometry
    with pytest.raises(AttributeError):
        machine.baxis
    with pytest.raises(AttributeError):
        machine.raxis
    with pytest.raises(AttributeError):
        machine.zaxis
    with pytest.raises(AttributeError):
        machine.rgeo
    with pytest.raises(AttributeError):
        machine.rlast
    with pytest.raises(RuntimeError):
        machine.quantity(1, "meter")
    machine.__repr__()
    machine.__str__()


def test_normalizations():
    lcfs = dex.LastClosedFluxSurface.Toroidal(0.05)
    machine = dex.Machine(
        geometry=dex.LarGeometry(2, 1.75, 0.5),
        qfactor=dex.UnityQfactor(lcfs),
        current=dex.LarCurrent(),
        bfield=dex.LarBfield(),
        species="Deuterium",
    )

    raxis_norm = machine.quantity(1, "NormMeter").to("meter").m
    assert isclose(raxis_norm, machine.raxis)
    baxis_norm = machine.quantity(1, "NormTesla").to("tesla").m
    assert isclose(baxis_norm, machine.baxis)


def test_gcmotion_units():
    geometry = dex.LarGeometry(3.5, 1.75, 0.5)
    LCFS = dex.LastClosedFluxSurface.Toroidal(geometry.psi_last)
    machine = dex.Machine(
        geometry=geometry,
        qfactor=dex.UnityQfactor(LCFS),
        current=dex.LarCurrent(),
        bfield=dex.LarBfield(),
        species="Deuterium",
    )

    Q = machine.quantity
    assert isclose(Q(1, "NormHertz").to("Hz").m, 167629580.2981652)
    assert isclose(Q(1e-5, "NormJoule").to("keV").m, 8.983897819104792)
    assert isclose(Q(1, "NormJoule").to("Joule").m, 1.4393791168013255e-10)
    assert isclose(Q(0.003, "NormHertz").to("kilohertz").m, 502.88874089449564)

    assert isclose(Q(geometry.baxis, "Tesla").to("NormTesla").m, 1)
    assert isclose(Q(1, "NormTesla").to("Tesla").m, geometry.baxis)
    assert isclose(Q(geometry.raxis, "Meter").to("NormMeter").m, 1)
    assert isclose(Q(1, "NormMeter").to("Meter").m, geometry.raxis)
