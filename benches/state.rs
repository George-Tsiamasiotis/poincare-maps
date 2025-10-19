use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};
use poincare_maps::*;

fn state_evaluation(c: &mut Criterion) {
    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1.0, 7.0, 0.0).unwrap(),
        Harmonic::from_dataset(&path, "akima", 1.0, 5.0, 0.0).unwrap(),
    ];
    let per = Perturbation::from_harmonics(harmonics);
    let psip_wall = qfactor.psip_wall;

    let initial = InitialConditions {
        t0: 0.0,
        theta0: 0.0,
        psip0: 0.0,
        rho0: 0.01,
        zeta0: 0.0,
        mu: 0.0,
    };

    let mut group = c.benchmark_group("State Evaluation");

    let mut state = State::new_init(&initial);
    let points = vec![(0.0, 0.5 * psip_wall)];
    let mut points_iter = points.iter().cycle();

    group.bench_function("fixed State::evaluate()", |b| {
        b.iter(|| {
            let next = points_iter.next().unwrap();
            state.theta = next.0;
            state.psip = next.1;
            state.evaluate(&qfactor, &current, &bfield, &per).unwrap();
        })
    });

    let mut state = State::new_init(&initial);
    let points = vec![(0.0, 0.1 * psip_wall), (3.14, 0.8 * psip_wall)];
    let mut points_iter = points.iter().cycle();

    group.bench_function("moving State::evaluate()", |b| {
        b.iter(|| {
            let next = points_iter.next().unwrap();
            state.theta = next.0;
            state.psip = next.1;
            state.evaluate(&qfactor, &current, &bfield, &per).unwrap();
        })
    });
}

criterion_group!(benches, state_evaluation);
criterion_main!(benches);
