import pytest
import dexter as dex
from math import isclose


def test_integrate_analytical(lar_machine_perturbed: dex.Machine):
    flux0 = dex.InitialFlux.Toroidal(0.01)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-4, 1e-6)
    particle = dex.Particle(initial)
    particle.integrate(lar_machine_perturbed, (0, 1e3))
    assert 10 < particle.steps_taken < 1000


def test_intersect_analytical(lar_machine_perturbed: dex.Machine):
    flux0 = dex.InitialFlux.Toroidal(0.025)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-6, 1e-7)
    particle = dex.Particle(initial)
    particle.intersect(lar_machine_perturbed, "ConstZeta", angle=0, turns=2)
    assert 10 < particle.steps_taken < 5000


def test_close_analytical(lar_machine: dex.Machine):
    flux0 = dex.InitialFlux.Toroidal(0.01)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-4, 1e-6)
    particle = dex.Particle(initial)
    particle.close(lar_machine)
    assert 10 < particle.steps_taken < 5000


def test_classify_analytical(lar_machine: dex.Machine):
    flux0 = dex.InitialFlux.Toroidal(0.01)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-4, 1e-6)
    particle = dex.Particle(initial)
    particle.classify(lar_machine)
