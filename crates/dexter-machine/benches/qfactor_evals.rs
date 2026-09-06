//! Benchmark for the `Qfactor` objects' evaluation methods.

#![allow(unused_results)]

use criterion::{Criterion, criterion_group, criterion_main};
use dexter_machine::extract::TEST_NETCDF_PATH;
use dexter_machine::*;

use std::hint::black_box;
use std::path::PathBuf;

fn qfactor_evals(c: &mut Criterion) {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let acc = &mut Accelerator::new();
    let psi = MagneticFlux::Toroidal(0.01);

    let unity = UnityQfactor::new(lcfs);
    let parabolic = ParabolicQfactor::new(1.1, 3.9, lcfs);
    let nc_qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Cubic)
        .build()
        .unwrap();

    // ===========================================================================================

    let mut group = c.benchmark_group("Qfactor q(ψ) evaluation");

    group.bench_with_input("Unity", &psi, |b, &psi| {
        b.iter(|| unity.eval_q(black_box(psi), black_box(acc)).unwrap());
    });
    group.bench_with_input("Parabolic", &psi, |b, &psi| {
        b.iter(|| parabolic.eval_q(black_box(psi), black_box(acc)).unwrap());
    });
    group.bench_with_input("NcQfactor", &psi, |b, &psi| {
        b.iter(|| nc_qfactor.eval_q(black_box(psi), black_box(acc)).unwrap());
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Qfactor ψ(ψp) evaluation");

    group.bench_with_input("Unity", &psi, |b, &psi| {
        b.iter(|| unity.eval_other(black_box(psi), black_box(acc)).unwrap());
    });
    group.bench_with_input("Parabolic", &psi, |b, &psi| {
        b.iter(|| {
            parabolic
                .eval_other(black_box(psi), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcQfactor", &psi, |b, &psi| {
        b.iter(|| {
            nc_qfactor
                .eval_other(black_box(psi), black_box(acc))
                .unwrap()
        });
    });
    group.finish();
}

criterion_group!(benches, qfactor_evals);
criterion_main!(benches);
