use std::path::PathBuf;

use equilibrium::*;

pub fn create_equilibrium() -> (Qfactor, Current, Bfield, Perturbation) {
    let path = PathBuf::from("../data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Current::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 7.0, 1.0, 0.0).unwrap(),
        Harmonic::from_dataset(&path, "akima", 8.0, 1.0, 0.0).unwrap(),
    ];
    let per = Perturbation::from_harmonics(&harmonics);

    (qfactor, current, bfield, per)
}
