use crate::*;
use pyo3::prelude::*;

use crate::error::DexterError;

pub fn resolve_interpolation_1d_type(interp_type: String) -> PyResult<Interpolation1dType> {
    use Interpolation1dType::*;
    Ok(match interp_type.to_lowercase().as_str() {
        "linear" => Linear,
        "cubic" => Cubic,
        "akima" => Akima,
        "cubicperiodic" => CubicPeriodic,
        "akimaperiodic" => AkimaPeriodic,
        "steffen" => Steffen,
        _ => {
            return Err(DexterError::InvalidInterpolationType(
                concat!(
                    "Supported 1D interpolation types are",
                    "'Linear', 'Cubic', 'CubicPeriodic', 'Akima', 'AkimaPeriodic' and 'Steffen'",
                )
                .into(),
            )
            .into());
        }
    })
}
