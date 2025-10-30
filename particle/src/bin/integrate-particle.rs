use std::{path::PathBuf, process::exit};

use equilibrium::{Bfield, Current, Harmonic, Perturbation, Qfactor};
use particle::*;

fn main() {
    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1.0, 7.0, 0.0).unwrap(),
        Harmonic::from_dataset(&path, "akima", 1.0, 9.0, 0.0).unwrap(),
    ];
    let per = Perturbation::from_harmonics(&harmonics);

    let psip_wall = qfactor.psip_wall;

    let initial = InitialConditions {
        t0: 0.0,
        theta0: 1.0,
        psip0: 0.5 * psip_wall,
        rho0: 0.005,
        zeta0: 0.0,
        mu: 0.0,
    };

    let mut particle = Particle::new(&initial);
    match particle.integrate(&qfactor, &bfield, &current, &per, (0.0, 100000.0)) {
        Ok(_) => dbg!(&particle),
        Err(err) => {
            eprintln!("{}", err);
            exit(2);
        }
    };
}
