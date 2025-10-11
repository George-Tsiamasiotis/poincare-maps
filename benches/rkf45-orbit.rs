use criterion::{Criterion, criterion_group, criterion_main};
use std::{path::PathBuf, time::Duration};

use poincare_maps::*;

fn rkf45_orbit(c: &mut Criterion) {
    #[cfg(feature = "rk45")]
    compile_error!("Feature rk45 must be disabled for this bench.");

    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let psip_wall = qfactor.psip_wall;

    let initial = InitialConditions::new(0.0, 0.0, 0.5 * psip_wall, 0.01, 0.0, 0.0);
    let mut particle = Particle::new(&initial);
    c.bench_function("rkf45 orbit", |b| {
        b.iter(|| particle.run_ode(&qfactor, &bfield, &current, (0.0, 20000.0), 0))
    });
    dbg!(&particle);
}

criterion_group!(benches, rkf45_orbit);
criterion_main!(benches);
