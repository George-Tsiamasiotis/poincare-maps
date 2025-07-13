extern crate gsl_splines;
extern crate ndarray;

use gsl_splines::{Accelerator, InterpolationType, Spline};
use ndarray::Array1;

fn main() {
    // Data creation
    let xdata = Array1::linspace(0.0, 3.0, 100);
    let ydata1 = xdata.sin();
    let ydata2 = xdata.cos();

    // Interpolation type and Accelerator
    let typ = InterpolationType::Cubic;
    let acc = Accelerator::new();

    // Spline creation
    let mut spline1 = Spline::build(typ, &xdata, &ydata1, acc.clone()).unwrap();
    let mut spline2 = Spline::build(typ, &xdata, &ydata2, acc.clone()).unwrap();

    // 2 Different splines evaluating on a different point.
    for x in Array1::linspace(1.0, 1.0001, 9) {
        spline1.eval(x).unwrap();
        spline2.eval(x + 1.0).unwrap();
    }

    println!("{:?}", spline1.acc);
    println!("{:?}", spline2.acc);

    spline1.reset_acc();
    spline2.reset_acc();
    println!();

    // 2 Different splines evaluating on the same point
    for x in Array1::linspace(1.0, 1.0001, 9) {
        spline1.eval(x).unwrap();
        spline2.eval(x).unwrap();
    }

    println!("{:?}", spline1.acc);
    println!("{:?}", spline2.acc);
}
