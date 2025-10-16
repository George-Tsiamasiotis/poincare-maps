use std::{f64::consts::PI, path::PathBuf};

use ndarray::Array1;
use poincare_maps::*;

fn setup_poincare() -> Result<(Bfield, Qfactor, Current, Perturbation, Poincare)> {
    let path = PathBuf::from("./data.nc");
    let bfield = Bfield::from_dataset(&path, "bicubic")?;
    let qfactor = Qfactor::from_dataset(&path, "akima")?;
    let current = Current::from_dataset(&path, "akima")?;
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1.0, 7.0, 0.0)?,
        Harmonic::from_dataset(&path, "akima", 1.0, 9.0, 0.0)?,
    ];
    let per = Perturbation::from_harmonics(harmonics);
    let psip_wall = qfactor.psip_wall;
    let n = 7;
    // Include some particles outside the wall.
    // Values are [0, 0.33, 0.66, 1, 1.33, 1.66, 2]. Only 0.33 and 0.66 should pass
    let psips = Array1::linspace(0.0, 2.0 * psip_wall, n);
    let mut map = Poincare::default();
    for psip in psips.iter() {
        let initial = InitialConditions::new(0.0, 0.0, *psip, 0.01, 0.0, 0.0);
        let particle = Particle::new(&initial);
        map.add_particle(&particle);
    }
    Ok((bfield, qfactor, current, per, map))
}

#[test]
fn test_poincare_theta() -> Result<()> {
    let (bfield, qfactor, current, per, mut map) = setup_poincare()?;

    let turns = 2;
    map.run(&qfactor, &bfield, &current, &per, "theta", PI, turns)?;

    let _ = map.get_particles();
    let _ = map.get_angles();
    let _ = map.get_fluxes();

    assert_eq!(map.completed_particles, 2);
    assert_eq!(map.escpaped_particles, 5);

    Ok(())
}

#[test]
fn test_poincare_zeta() -> Result<()> {
    let (bfield, qfactor, current, per, mut map) = setup_poincare()?;

    let turns = 2;
    map.run(&qfactor, &bfield, &current, &per, "zeta", PI, turns)?;

    let _ = map.get_particles();
    let _ = map.get_angles();
    let _ = map.get_fluxes();

    assert_eq!(map.completed_particles, 2);
    assert_eq!(map.escpaped_particles, 5);

    Ok(())
}
