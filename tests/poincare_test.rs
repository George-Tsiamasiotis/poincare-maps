use std::{f64::consts::PI, path::PathBuf};

use ndarray::Array1;
use poincare_maps::*;

fn setup_poincare(
    params: PoincareParameters,
) -> Result<(Bfield, Qfactor, Current, Perturbation, Poincare)> {
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
    // Values are [0, 0.33, 0.66, 1, 1.33, 1.66, 2]. Only 0, 0.33 and 0.66 should pass
    let psips = Array1::linspace(0.0, 2.0 * psip_wall, n);
    let mut map = Poincare::new(params);
    for psip in psips.iter() {
        let particle = Particle::new(0.0, 0.0, *psip, 0.01, 0.0, 0.0);
        map.add_particle(&particle);
    }
    Ok((bfield, qfactor, current, per, map))
}

#[test]
fn test_poincare_theta() -> Result<()> {
    let params = PoincareParameters::new("theta", PI, 2);
    let (bfield, qfactor, current, per, mut map) = setup_poincare(params)?;

    map.run(&qfactor, &bfield, &current, &per)?;

    let _ = map.get_particles();

    assert_eq!(map.completed_particles, 3);
    assert_eq!(map.escpaped_particles, 4);

    Ok(())
}

#[test]
fn test_poincare_zeta() -> Result<()> {
    let params = PoincareParameters::new("zeta", PI, 2);
    let (bfield, qfactor, current, per, mut map) = setup_poincare(params)?;

    map.run(&qfactor, &bfield, &current, &per)?;

    let _ = map.get_particles();

    assert_eq!(map.completed_particles, 3);
    assert_eq!(map.escpaped_particles, 4);

    Ok(())
}
