use std::f64::consts::PI;
use std::path::PathBuf;
use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use equilibrium::*;
use particle::State;
use particle::*;

fn state_evaluation(c: &mut Criterion) {
    let path = PathBuf::from("../data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Currents::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1, 7).unwrap(),
        Harmonic::from_dataset(&path, "akima", 1, 9).unwrap(),
    ];
    let per = Perturbation::from_harmonics(&harmonics);
    let psip_wall = qfactor.psip_wall();

    let initial = InitialConditions {
        time0: 0.0,
        theta0: 0.0,
        psip0: 0.0,
        rho0: 0.01,
        zeta0: 0.0,
        mu: 0.0,
    };

    let mut state = State::from_initial(&initial);
    let points = vec![(10.0 * PI, 0.1 * psip_wall), (15.0 * PI, 0.8 * psip_wall)];
    let mut points_iter = points.iter().cycle();

    let mut group = c.benchmark_group("State Evaluation");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(3));

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
