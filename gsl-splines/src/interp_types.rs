use crate::{RgslInterp2dType, RgslInterpType};

/// Spline types supported by GSL.
///
/// This type is a thin wrapper around `rgsl`'s `InterpType`.
///
/// ## Example
/// ```
/// # use gsl_splines::InterpolationType;
/// #
/// # fn main() {
/// let spline_type = gsl_splines::InterpolationType::Akima;
/// # }
/// ```
/// Descriptions are taken directly from GSL's
/// [interpolation](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_type) page
#[derive(Debug, Clone, Copy)]
pub enum InterpolationType {
    /// Linear interpolation. This interpolation method does not require any additional memory.
    Linear,

    /// Polynomial interpolation. This method should only be used for interpolating small number
    /// of points because polynomial interpolation introduces large oscillations, even for
    /// well-behaved datasets. The number of terms in the interpolating polynomial is equal to
    /// the number of points.
    Polynomial,

    /// Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on
    /// each interval, with matching first and second derivatives at the supplied data-points.
    /// The second derivative is chosen to be zero at the first point and last point.
    Cubic,

    /// Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic
    /// on each interval, with matching first and second derivatives at the supplied data-points.
    /// The derivatives at the first and last points are also matched. Note that the last point
    /// in the data must have the same y-value as the first point, otherwise the resulting
    /// periodic interpolation will have a discontinuity at the boundary.
    CubicPeriodic,

    /// Non-rounded Akima spline with natural boundary conditions. This method uses the
    /// non-rounded corner algorithm of Wodicka.
    Akima,

    /// Non-rounded Akima spline with periodic boundary conditions. This method uses the
    /// non-rounded corner algorithm of Wodicka.
    AkimaPeriodic,
    ///Steffen’s method guarantees the monotonicity of the interpolating function between the
    ///given data points. Therefore, minima and maxima can only occur exactly at the data points,
    ///and there can never be spurious oscillations between data points. The interpolation
    /// function is piecewise cubic in each interval. The resulting curve and its first derivative
    /// are guaranteed to be continuous, but the second derivative may be discontinuous.
    Steffen,
}

impl From<InterpolationType> for RgslInterpType {
    /// Get the corresponding gsl_interp_type.
    fn from(value: InterpolationType) -> RgslInterpType {
        use InterpolationType::*;

        match value {
            Linear => RgslInterpType::linear(),
            Polynomial => RgslInterpType::polynomial(),
            Cubic => RgslInterpType::cspline(),
            CubicPeriodic => RgslInterpType::cspline_periodic(),
            Akima => RgslInterpType::akima(),
            AkimaPeriodic => RgslInterpType::akima_periodic(),
            Steffen => unimplemented!("Steffen splines are implemented in GSL, but not rgsl yet."),
        }
    }
}

/// Spline2D types supported by GSL.
///
/// This type is a thin wrapper around `rgsl`'s `Interp2dType`.
///
/// ## Example
/// ```
/// # use gsl_splines::Interpolation2dType;
/// #
/// # fn main() {
/// let spline2d_type = gsl_splines::Interpolation2dType::Bicubic;
/// # }
/// ```
/// Descriptions are taken directly from GSL's
/// [interpolation](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp2d_type) page
#[derive(Debug, Clone, Copy)]
pub enum Interpolation2dType {
    /// Bilinear interpolation. This interpolation method does not require any additional memory.
    Bilinear,
    /// Bicubic interpolation.
    Bicubic,
}

impl From<Interpolation2dType> for RgslInterp2dType {
    /// Get the corresponding gsl_interp2d_type.
    fn from(value: Interpolation2dType) -> RgslInterp2dType {
        use Interpolation2dType::*;

        match value {
            Bilinear => RgslInterp2dType::bilinear(),
            Bicubic => RgslInterp2dType::bicubic(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_interpolationtype_from() {
        let _: RgslInterpType = InterpolationType::Linear.into();
        let _: RgslInterpType = InterpolationType::Polynomial.into();
        let _: RgslInterpType = InterpolationType::Cubic.into();
        let _: RgslInterpType = InterpolationType::CubicPeriodic.into();
        let _: RgslInterpType = InterpolationType::Akima.into();
        let _: RgslInterpType = InterpolationType::AkimaPeriodic.into();
    }

    #[test]
    fn test_interpolation2dtype_from() {
        let _: RgslInterp2dType = Interpolation2dType::Bilinear.into();
        let _: RgslInterp2dType = Interpolation2dType::Bicubic.into();
    }
}
