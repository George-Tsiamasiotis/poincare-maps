//! Wrappers around `rgsl`'s `eval_<>_e()' functions. These functions are prefered over their 'eval_<>()'
//! alternatives, since they can catch GSL's errors. The overhead seems to be negligible.

use std::f64;

use crate::Result;
use crate::{SplineError, spline::Spline};

impl Spline {
    #[inline]
    /// Returns the interpolated value of `y` for a given point `x`.
    ///
    /// ## Errors
    ///
    /// Returns an error when `x` is outside the range of the supplied data points.
    pub fn eval(&mut self, x: f64) -> Result<f64> {
        match self
            .gsl_spline
            .eval_e(x, &mut self.acc.borrow_mut().gsl_iterp_accel)
        {
            Ok(val) => Ok(val),
            Err(err) => Err(SplineError::GSLInputDomainError { err }),
        }
    }

    #[inline]
    /// Returns the derivative of an interpolated function for a given point `x`.
    ///
    /// ## Errors
    ///
    /// Returns an error when `x` is outside the range of the supplied data points.
    pub fn eval_deriv(&mut self, x: f64) -> Result<f64> {
        match self
            .gsl_spline
            .eval_deriv_e(x, &mut self.acc.borrow_mut().gsl_iterp_accel)
        {
            Ok(val) => Ok(val),
            Err(err) => Err(SplineError::GSLInputDomainError { err }),
        }
    }

    #[inline]
    /// Returns the second derivative of an interpolated function for a given point `x`.
    ///
    /// ## Errors
    ///
    /// Returns an error when `x` is outside the range of the supplied data points.
    pub fn eval_deriv2(&mut self, x: f64) -> Result<f64> {
        match self
            .gsl_spline
            .eval_deriv2_e(x, &mut self.acc.borrow_mut().gsl_iterp_accel)
        {
            Ok(val) => Ok(val),
            Err(err) => Err(SplineError::GSLInputDomainError { err }),
        }
    }

    #[inline]
    /// Returns the numerical integral of an interpolated function over the range [a, b].
    ///
    /// ## Errors
    /// Returns an Error if `a>b`, or if either `a` or `b` are outside the range of the supplied data points.
    pub fn eval_integ(&mut self, a: f64, b: f64) -> Result<f64> {
        match self
            .gsl_spline
            .eval_integ_e(a, b, &mut self.acc.borrow_mut().gsl_iterp_accel)
        {
            Ok(val) => Ok(val),
            Err(err) => Err(SplineError::GSLInputDomainError { err }),
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{Accelerator, InterpolationType, RgslValue, Spline, SplineError};
    use ndarray::Array1;

    fn some_spline() -> Spline {
        let xdata = Array1::linspace(0.0, 3.0, 10);
        let ydata = xdata.pow2();
        let typ = InterpolationType::CubicPeriodic;
        let acc = Accelerator::new();
        Spline::build(typ, &xdata, &ydata, acc).unwrap()
    }

    #[test]
    fn test_eval() {
        assert!(some_spline().eval(2.0).is_ok());
        assert!(matches!(
            some_spline().eval(200.0).unwrap_err(),
            SplineError::GSLInputDomainError { .. }
        ))
    }

    #[test]
    fn test_deriv() {
        assert!(some_spline().eval_deriv(2.0).is_ok());
        assert!(matches!(
            some_spline().eval_deriv(200.0).unwrap_err(),
            SplineError::GSLInputDomainError { .. }
        ))
    }

    #[test]
    fn test_deriv2() {
        assert!(some_spline().eval_deriv2(2.0).is_ok());
        assert!(matches!(
            some_spline().eval_deriv2(200.0).unwrap_err(),
            SplineError::GSLInputDomainError { .. }
        ))
    }

    #[test]
    fn test_integ() {
        assert!(some_spline().eval_integ(1.0, 2.0).is_ok());
        // All checks that GSL does
        assert!(matches!(
            some_spline().eval_integ(2.0, 1.0).unwrap_err(),
            SplineError::GSLInputDomainError { .. }
        ));
        assert!(matches!(
            some_spline().eval_integ(-2.0, 1.0).unwrap_err(),
            SplineError::GSLInputDomainError { .. }
        ));
        assert!(matches!(
            some_spline().eval_integ(2.0, 100.0).unwrap_err(),
            SplineError::GSLInputDomainError { .. }
        ));
        // Check that it prints
        let _ = format!(
            "{:?}",
            SplineError::GSLInputDomainError {
                err: RgslValue::Domain
            }
        );
    }
}
