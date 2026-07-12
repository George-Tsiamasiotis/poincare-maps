//! Benchmark for the equilibrium objects evaluation methods.

#![allow(unused_results)]

use criterion::{Criterion, criterion_group, criterion_main};
use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use rsl_interpolation::{Accelerator, Cache};

use std::path::PathBuf;

fn evaluations_benchmark(c: &mut Criterion) {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let mut psi_acc = Accelerator::new();
    let mut theta_acc = Accelerator::new();
    let mut cache = Cache::new();
    let (t, psi, theta, zeta) = (10000.0, 0.01, 1.0, 4.0);

    // ===========================================================================================

    let unity_qfactor = UnityQfactor::new(lcfs);
    let parabolic_qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
    let nc_qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();

    let mut group = c.benchmark_group("Qfactor q(ψ) evaluation");

    group.bench_with_input("UnityQfactor", &psi, |b, &psi| {
        b.iter(|| unity_qfactor.q_of_psi(psi, &mut psi_acc));
    });
    group.bench_with_input("ParabolicQfactor", &psi, |b, &psi| {
        b.iter(|| parabolic_qfactor.q_of_psi(psi, &mut psi_acc));
    });
    group.bench_with_input("NcQfactor", &psi, |b, &psi| {
        b.iter(|| nc_qfactor.q_of_psi(psi, &mut psi_acc));
    });
    group.finish();

    // ===========================================================================================

    let lar_current = LarCurrent::new();
    let nc_current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();

    let mut group = c.benchmark_group("Current g(ψ) evaluation");

    group.bench_with_input("LarCurrent", &psi, |b, &psi| {
        b.iter(|| lar_current.g_of_psi(psi, &mut psi_acc));
    });
    group.bench_with_input("NcCurrent", &psi, |b, &psi| {
        b.iter(|| nc_current.g_of_psi(psi, &mut psi_acc));
    });
    group.finish();

    // ===========================================================================================

    let lar_bfield = LarBfield::new();
    let nc_bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();

    let mut group = c.benchmark_group("Bfield B(ψ, θ) evaluation");

    group.bench_with_input("LarBfield", &(psi, theta), |b, &(psi, theta)| {
        b.iter(|| lar_bfield.b_of_psi(psi, theta, &mut psi_acc, &mut theta_acc, &mut cache));
    });
    group.bench_with_input("NcBfield", &(psi, theta), |b, &(psi, theta)| {
        b.iter(|| nc_bfield.b_of_psi(psi, theta, &mut psi_acc, &mut theta_acc, &mut cache));
    });
    group.finish();

    // ===========================================================================================

    let flute_mode = FluteMode::new(1e-3, lcfs, 1, 2, 0.0);
    let nc_flute_mode = NcFluteModeBuilder::new(&path, "steffen", 2, 1)
        .with_phase_method(PhaseMethod::Interpolation)
        .build()
        .unwrap();
    let mut flute_mode_cache = flute_mode.generate_cache();
    let mut nc_flute_mode_cache = nc_flute_mode.generate_cache();

    let mut group = c.benchmark_group("Flute mode α(ψ, θ , ζ, t) evaluation");

    group.bench_with_input(
        "Analytical flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| flute_mode.h_of_psi(psi, theta, zeta, t, &mut flute_mode_cache));
        },
    );
    group.bench_with_input(
        "Nc flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| nc_flute_mode.h_of_psi(psi, theta, zeta, t, &mut nc_flute_mode_cache));
        },
    );
    group.finish();
}

criterion_group!(benches, evaluations_benchmark);
criterion_main!(benches);
