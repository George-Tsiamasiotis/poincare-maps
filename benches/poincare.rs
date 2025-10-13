use std::time::Duration;
use std::{f64::consts::PI, path::PathBuf};

use criterion::Criterion;
use criterion::{criterion_group, criterion_main};
use ndarray::Array1;
use poincare_maps::*;

fn poincare(c: &mut Criterion) {
    #[cfg(feature = "rk45")]
    compile_error!("Feature rk45 must be disabled for this bench.");

    let path = PathBuf::from("./data.nc");
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();
    let per = Perturbation::from_dataset(&path, "akima", 1.0, -8.0).unwrap();
    let psip_wall = qfactor.psip_wall;

    let n = 20;
    let turns = 100;
    let psips = Array1::linspace(0.05 * psip_wall, 0.95 * psip_wall, n);

    let mut map = Poincare::new();

    for psip in psips.iter() {
        let initial = InitialConditions::new(0.0, 0.0, *psip, 0.01, 0.0, 0.0);
        let particle = Particle::new(&initial);
        map.add_particle(&particle);
    }

    let mut g = c.benchmark_group("Poincare Map Calculation");
    g.measurement_time(Duration::from_secs(30));
    g.bench_function("poincare", |b| {
        b.iter(|| map.run(&qfactor, &bfield, &current, &per, "theta", PI, turns))
    });
}

criterion_group!(benches, poincare);
criterion_main!(benches);
