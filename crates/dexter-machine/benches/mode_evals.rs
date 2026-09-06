//! Benchmark for the `Mode` objects' evaluation methods.

#![allow(unused_results)]

use criterion::{Criterion, criterion_group, criterion_main};
use dexter_machine::extract::TEST_NETCDF_PATH;
use dexter_machine::*;

use std::hint::black_box;
use std::path::PathBuf;

fn mode_evals(c: &mut Criterion) {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let (t, psi, theta, zeta) = (10000.0, MagneticFlux::Toroidal(0.01), 1.0, 4.0);

    let flute_mode = FluteMode::new(1e-3, lcfs, 1, 2, 0.0);
    let nc_flute_mode = NcFluteModeBuilder::new(&path, Interpolation1dType::Steffen, 2, 1)
        .with_phase_method(PhaseMethod::Interpolation)
        .build()
        .unwrap();
    let mut flute_mode_cache = flute_mode.generate_cache();
    let mut nc_flute_mode_cache = nc_flute_mode.generate_cache();

    // ===========================================================================================

    let mut group = c.benchmark_group("Flute mode m(ψ, θ , ζ, t) evaluation");

    group.bench_with_input(
        "Analytical flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                flute_mode
                    .eval_m(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.bench_with_input(
        "Nc flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                nc_flute_mode
                    .eval_m(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut nc_flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Flute mode dm(ψ, θ , ζ, t)/dψ evaluation");

    group.bench_with_input(
        "Analytical flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                flute_mode
                    .eval_deriv_flux(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.bench_with_input(
        "Nc flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                nc_flute_mode
                    .eval_deriv_flux(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut nc_flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Flute mode dm(ψ, θ , ζ, t)/dθ evaluation");

    group.bench_with_input(
        "Analytical flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                flute_mode
                    .eval_deriv_theta(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.bench_with_input(
        "Nc flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                nc_flute_mode
                    .eval_deriv_theta(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut nc_flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Flute mode dm(ψ, θ , ζ, t)/dζ evaluation");

    group.bench_with_input(
        "Analytical flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                flute_mode
                    .eval_deriv_zeta(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.bench_with_input(
        "Nc flute mode",
        &(psi, theta, zeta, t),
        |b, &(psi, theta, zeta, t)| {
            b.iter(|| {
                nc_flute_mode
                    .eval_deriv_zeta(
                        black_box(psi),
                        black_box(theta),
                        black_box(zeta),
                        black_box(t),
                        black_box(&mut nc_flute_mode_cache),
                    )
                    .unwrap()
            });
        },
    );
    group.finish();
}

criterion_group!(benches, mode_evals);
criterion_main!(benches);
