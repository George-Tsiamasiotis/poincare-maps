use std::f64::consts::PI;
use std::path::PathBuf;
use std::time::Duration;

use crate::PoincareSection::ConstTheta;
use criterion::{Criterion, criterion_group, criterion_main};
use equilibrium::*;
use particle::*;

fn state_evaluation(c: &mut Criterion) {
    let path = PathBuf::from("../data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let currents = Currents::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1, 7).unwrap(),
        Harmonic::from_dataset(&path, "akima", 1, 9).unwrap(),
    ];
    let perturbation = Perturbation::from_harmonics(&harmonics);

    let psip_wall = qfactor.psip_wall();

    let initial = InitialConditions {
        time0: 0.0,
        theta0: 1.0,
        psip0: 0.3 * psip_wall,
        rho0: 0.005,
        zeta0: 0.0,
        mu: 0.0,
    };

    let mut particle = Particle::new(&initial);
    let params = MappingParameters::new(ConstTheta, PI, 50);

    let mut group = c.benchmark_group("Particle Map");
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(5));

    group.bench_function("Particle Map", |b| {
        b.iter(|| particle.map(&qfactor, &bfield, &currents, &perturbation, &params))
    });
}

criterion_group!(benches, state_evaluation);
criterion_main!(benches);
