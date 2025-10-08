use poincare_maps::*;
use std::path::PathBuf;

fn main() {
    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();

    // Normal passing particle
    let initial = InitialConditions {
        t0: 0.0,
        theta0: 3.14,
        psip0: 0.05,
        rho0: 0.05,
        zeta0: 0.1,
        mu: 1e-4,
    };
    let mut particle = Particle::new(initial);
    particle
        .run_ode(&qfactor, &bfield, &current, (0.0, 200.0), 5000)
        .unwrap();
    dbg!(&particle);
}
