mod common;

use std::f64::consts::PI;

use ndarray::Array1;
use particle::{MappingParameters, PoincareSection};
use poincare::*;

use crate::common::create_equilibrium;

#[test]
fn test_normal_particle_int() {
    let (qfactor, current, bfield, per) = create_equilibrium();
    let psip_wall = qfactor.psip_wall();

    let num = 10;
    let init = PoincareInit::build(
        Array1::from_elem(num, 0.0).as_slice().unwrap(),
        Array1::linspace(0.0, psip_wall * 0.8, num)
            .as_slice()
            .unwrap(),
        Array1::from_elem(num, 0.001).as_slice().unwrap(),
        Array1::from_elem(num, 0.0).as_slice().unwrap(),
        Array1::from_elem(num, 0.0).as_slice().unwrap(),
    )
    .unwrap();

    let intersections = 10;
    let params = MappingParameters::new(PoincareSection::ConstTheta, PI, intersections);

    let mut p = Poincare::new(init, params);
    p.run(&qfactor, &bfield, &current, &per).unwrap();

    assert_eq!(p.angles().shape(), &[num, intersections]);
    assert_eq!(p.fluxes().shape(), &[num, intersections]);
}
