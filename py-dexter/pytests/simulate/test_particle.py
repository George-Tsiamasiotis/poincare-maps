import pytest
import numpy as np
import dexter as dex
from math import isclose, isfinite


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


def test_getters(lar_machine: dex.Machine):
    flux0 = dex.InitialFlux.Toroidal(0.01)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-4, 1e-6)
    particle = dex.Particle(initial)
    particle.close(lar_machine)
    particle.classify(lar_machine)
    assert 10 < particle.steps_taken < 5000

    assert isinstance(particle.initial_conditions, dex.InitialConditions)
    assert particle.integration_status == "ClosedPeriods(1)"
    assert particle.steps_taken == particle.steps_stored - 1
    assert isinstance(particle.duration, str)
    assert particle.initial_energy is not None and isfinite(particle.initial_energy)
    assert particle.final_energy is not None and isfinite(particle.final_energy)
    assert particle.energy_var is not None and isfinite(particle.energy_var)
    assert particle.energy_pzeta_position == "Iota"
    assert particle.orbit_type == "TrappedConfined"
    assert isfinite(particle.omega_theta)
    assert isfinite(particle.omega_zeta)
    assert isfinite(particle.qkinetic)
    particle.print_caches()
    assert particle.flux_cache_hits == 0
    assert particle.flux_cache_misses == 0
    assert particle.theta_cache_hits == 0
    assert particle.theta_cache_misses == 0
    assert particle.mode_cache_hits == 0
    assert particle.mode_cache_misses == 0

    assert len(particle.t_array) == particle.steps_stored
    assert len(particle.psi_array) == particle.steps_stored
    assert len(particle.psip_array) == particle.steps_stored
    assert len(particle.theta_array) == particle.steps_stored
    assert len(particle.zeta_array) == particle.steps_stored
    assert len(particle.rho_array) == particle.steps_stored
    assert len(particle.mu_array) == particle.steps_stored
    assert len(particle.ptheta_array) == particle.steps_stored
    assert len(particle.pzeta_array) == particle.steps_stored
    assert len(particle.energy_array) == particle.steps_stored
    particle.discard_arrays()
    assert len(particle.t_array) == 0
    particle.__str__()
    particle.__repr__()
