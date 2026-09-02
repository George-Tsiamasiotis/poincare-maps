//! Provides functions for instantiating `Box<dyn Interpolation>` objects.

use rsl_interpolation::{
    AkimaInterpolator, AkimaPeriodicInterpolator, BicubicInterpolator, BilinearInterpolator,
    BuildInterpolator, BuildInterpolator2d, CubicInterpolator, CubicPeriodicInterpolator,
    Interpolation, Interpolation2d, InterpolationError, LinearInterpolator, SteffenInterpolator,
};

/// Type alias for thread-safe and cloneable [`Interpolation`] objects.
pub type DynInterpolator = Box<dyn Interpolation>;

/// Type alias for thread-safe and cloneable [`Interpolation2d`] objects.
pub type DynInterpolator2d = Box<dyn Interpolation2d>;

/// Available 1D interpolation types.
///
/// Types are provided from the [`rsl_interpolation`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation1dType {
    /// Linear interpolation.
    ///
    /// This interpolation method does not require any additional memory.
    Linear,
    /// Cubic spline with natural boundary conditions.
    ///
    /// The resulting curve is piecewise cubic on each interval, with matching first and second
    /// derivatives at the supplied data-points. The second derivative is chosen to be zero at the
    /// first point and last point.
    Cubic,
    /// Cubic spline with periodic boundary conditions.
    ///
    /// The resulting curve is piecewise cubic on each interval, with matching first and second
    /// derivatives at the supplied data-points. The derivatives at the first and last points are
    /// also matched. Note that the last point in the data must have the same y-value as the first
    /// point, otherwise the resulting periodic interpolation will have a discontinuity at the
    /// boundary.
    CubicPeriodic,
    /// Non-rounded Akima spline with natural boundary conditions.
    ///
    /// This method uses the non-rounded corner algorithm of Wodicka.
    Akima,
    /// Non-rounded Akima spline with periodic boundary conditions.
    ///
    /// This method uses the non-rounded corner algorithm of Wodicka.
    AkimaPeriodic,
    /// Steffen interpolation.
    ///
    /// Steffen’s method guarantees the monotonicity of the interpolating function between the
    /// given data points. Therefore, minima and maxima can only occur exactly at the data
    /// points, and there can never be spurious oscillations between data points. The interpolated
    /// function is piecewise cubic in each interval. The resulting curve and its first derivative
    /// are guaranteed to be continuous, but the second derivative may be discontinuous.
    Steffen,
}

/// Available 2D interpolation types.
///
/// Types are provided from the [`rsl_interpolation`].
///
/// ## References
///
/// Numerical Algorithms with C - Gisela Engeln-Mullges, Frank Uhlig - 1996 -
/// Algorithm 10.1, pg 254.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation2dType {
    /// Bilinear interpolation.
    ///
    /// This interpolation method does not require any additional memory.
    Bilinear,
    /// Bicubic Interpolation.
    Bicubic,
}

/// Constructs a [`DynInterpolator`] of interpolation type `typ` and data arrays `xa` and `ya`.
///
/// # Errors
///
/// Returns an [`InterpolationError`] if the interpolator's `build` method fails. See
/// interpolator's documentation for possible errors.
pub fn make_interp(
    typ: Interpolation1dType,
    xa: &[f64],
    ya: &[f64],
) -> Result<DynInterpolator, InterpolationError> {
    use Interpolation1dType::*;
    let interp: DynInterpolator = match typ {
        Linear => Box::new(LinearInterpolator::build(xa, ya)?),
        Cubic => Box::new(CubicInterpolator::build(xa, ya)?),
        CubicPeriodic => Box::new(CubicPeriodicInterpolator::build(xa, ya)?),
        Akima => Box::new(AkimaInterpolator::build(xa, ya)?),
        AkimaPeriodic => Box::new(AkimaPeriodicInterpolator::build(xa, ya)?),
        Steffen => Box::new(SteffenInterpolator::build(xa, ya)?),
    };
    Ok(interp)
}

/// Constructs a [`DynInterpolator2d`] of interpolation type `typ` and data arrays `xa`, `ya` and `za`.
///
/// # Errors
///
/// Returns an [`InterpolationError`] if the interpolator's `build` method fails. See
/// interpolator's documentation for possible errors.
pub fn make_interp2d(
    typ: Interpolation2dType,
    xa: &[f64],
    ya: &[f64],
    za: &[f64],
) -> Result<DynInterpolator2d, InterpolationError> {
    use Interpolation2dType::*;
    let interp: DynInterpolator2d = match typ {
        Bilinear => Box::new(BilinearInterpolator::build(xa, ya, za)?),
        Bicubic => Box::new(BicubicInterpolator::build(xa, ya, za)?),
    };
    Ok(interp)
}
