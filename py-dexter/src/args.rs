use pyo3::types::PyTuple;

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

pub fn resolve_phase_method<'py>(arg: Bound<'py, PyAny>) -> Result<PhaseMethod> {
    use PhaseMethod::*;

    match arg.to_string().to_lowercase().as_str() {
        "zero" => return Ok(Zero),
        "average" => return Ok(Average),
        "interpolation" => return Ok(Interpolation),
        _ => (),
    }

    let tuple = match arg.cast::<PyTuple>() {
        Ok(tuple) => tuple,
        Err(_) => return Err(DexterError::InvalidPhaseMethod),
    };
    let string = tuple.get_item(0)?.extract::<String>()?.to_lowercase();
    let value = tuple.get_item(1)?.extract::<f64>()?;
    match string.as_str() {
        "custom" if value.is_finite() => Ok(Custom(value)),
        _ => Err(DexterError::InvalidPhaseMethod),
    }
}
