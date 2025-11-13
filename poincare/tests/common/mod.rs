use std::path::PathBuf;

use equilibrium::*;

pub fn create_equilibrium() -> (Qfactor, Currents, Bfield, Perturbation) {
    let path = PathBuf::from("../data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let current = Currents::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1, 7).unwrap(),
        Harmonic::from_dataset(&path, "akima", 1, 9).unwrap(),
    ];
    let per = Perturbation::from_harmonics(&harmonics);

    (qfactor, current, bfield, per)
}
