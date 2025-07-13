use std::f64;

use crate::Result;
use crate::{SplineError, spline::Spline};

/// Wrappers around `rgsl`'s `eval_<>_e()' functions. These functions are prefered over their
/// 'eval_<>()' alternatives, since they can catch GSL's errors. The overhead seems to be
/// negligible.
impl Spline {
    #[inline]
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
