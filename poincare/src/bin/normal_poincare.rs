use std::f64::consts::PI;
use std::process::exit;

use equilibrium::*;
use ndarray::Array1;
use particle::{MappingParameters, PoincareSection::*};
use poincare::*;

fn main() {
    let path = std::path::PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let currents = Currents::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 7.0, 1.0, 0.0).unwrap(),
        Harmonic::from_dataset(&path, "akima", 8.0, 1.0, 0.0).unwrap(),
    ];
    let perturbation = Perturbation::from_harmonics(&harmonics);
    let psip_wall = qfactor.psip_wall();

    let num = 100;
    let init = PoincareInit::build(
        Array1::from_elem(num, 0.0).as_slice().unwrap(),
        Array1::linspace(0.0, psip_wall, num).as_slice().unwrap(),
        Array1::from_elem(num, 0.01).as_slice().unwrap(),
        Array1::from_elem(num, 0.0).as_slice().unwrap(),
        Array1::from_elem(num, 0.0).as_slice().unwrap(),
    )
    .unwrap();

    let params = MappingParameters::new(ConstTheta, PI, 500);

    let mut p = Poincare::new(init, params);
    match p.run(&qfactor, &bfield, &currents, &perturbation) {
        Ok(_) => dbg!(&p),
        Err(err) => {
            eprintln!("{}", err);
            exit(2);
        }
    };
}
