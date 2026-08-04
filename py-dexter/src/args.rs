use crate::error::DexterError;
use crate::*;

pub fn resolve_interpolation_1d_type(interp_type: String) -> Result<Interpolation1dType> {
    use Interpolation1dType::*;
    Ok(match interp_type.to_lowercase().as_str() {
        "linear" => Linear,
        "cubic" => Cubic,
        "akima" => Akima,
        "cubicperiodic" => CubicPeriodic,
        "akimaperiodic" => AkimaPeriodic,
        "steffen" => Steffen,
        _ => return Err(DexterError::InvalidInterpolation1dType),
    })
}

pub fn resolve_interpolation_2d_type(interp_type: String) -> Result<Interpolation2dType> {
    use Interpolation2dType::*;
    Ok(match interp_type.to_lowercase().as_str() {
        "bilinear" => Bilinear,
        "bicubic" => Bicubic,
        _ => return Err(DexterError::InvalidInterpolation2dType),
    })
}
