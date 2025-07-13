extern crate gsl_splines;
extern crate ndarray;

use gsl_splines::{Accelerator, InterpolationType, Spline};
use ndarray::Array1;

fn main() {
    // Data creation
    let xdata = Array1::linspace(0.0, 3.0, 100);
    let ydata = xdata.pow2();

    // Interpolation type and Accelerator
    let typ = InterpolationType::Akima;
    let acc = Accelerator::new();

    // Spline creation
    let mut spline = Spline::build(typ, &xdata, &ydata, acc).unwrap();

    let _: f64 = spline.eval(1.0).unwrap();
    let _: f64 = spline.eval(1.3).unwrap();

    // Accelerator kicks in
    let _: f64 = spline.eval(2.0).unwrap();
    let _: f64 = spline.eval(2.0f64.next_up()).unwrap();

    println!("{spline:#?}");

    assert!(spline.eval(1000.0).is_err());
}
