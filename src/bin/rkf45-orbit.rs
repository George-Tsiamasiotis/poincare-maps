use std::path::PathBuf;

use poincare_maps::*;

fn main() {
    #[cfg(feature = "rk45")]
    compile_error!("Feature rk45 must be disabled for this bench.");

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

    let initial = InitialConditions::new(0.0, 0.0, 0.5 * psip_wall, 0.01, 0.0, 0.0);
    let mut particle = Particle::new(&initial);
    particle
        .run_ode(&qfactor, &bfield, &current, &per, (0.0, 20000.0), 0)
        .unwrap();

    dbg!(&particle);
}
