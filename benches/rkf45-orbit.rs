use criterion::{Criterion, criterion_group, criterion_main};
use std::{path::PathBuf, time::Duration};

use poincare_maps::*;

fn rkf45_orbit(c: &mut Criterion) {
    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1.0, 7.0, 0.0).unwrap(),
        // Harmonic::from_dataset(&path, "akima", 1.0, 5.0, 0.0)
    ];
    let per = Perturbation::from_harmonics(harmonics);
    let psip_wall = qfactor.psip_wall;

    let mut group = c.benchmark_group("rkf45 orbit");
    group.measurement_time(Duration::from_secs(10));

    group.bench_function("rkf45 orbit run()", |b| {
        b.iter(|| {
            let mut particle = Particle::new(0.0, 0.0, 0.5 * psip_wall, 0.01, 0.0, 0.0);
            particle.run_ode(&qfactor, &bfield, &current, &per, (0.0, 300.0))
        })
    });
}

criterion_group!(benches, rkf45_orbit);
criterion_main!(benches);
