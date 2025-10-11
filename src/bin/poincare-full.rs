use std::{f64::consts::PI, path::PathBuf};

use ndarray::Array1;
use poincare_maps::*;

fn main() -> Result<()> {
    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima")?;
    let current = Current::from_dataset(&path, "akima")?;
    let bfield = Bfield::from_dataset(&path, "bicubic")?;
    let psip_wall = qfactor.psip_wall;

    let n = 60;
    let turns = 200;
    let psips = Array1::linspace(0.05 * psip_wall, 0.95 * psip_wall, n);

    let mut map = Poincare::new();

    for psip in psips.iter() {
        let initial = InitialConditions::new(0.0, 0.0, *psip, 0.01, 0.0, 0.0);
        let particle = Particle::new(&initial);
        map.add_particle(&particle);
    }

    map.run(&qfactor, &bfield, &current, "theta", PI, turns)
        .unwrap();

    Ok(())
}
