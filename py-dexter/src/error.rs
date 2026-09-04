use dexter::dexter_machine::*;
use pyo3::{CastError, exceptions::PyException, prelude::*};

#[derive(Debug)]
pub enum DexterError {
    PyErr(String),
    CastError(String),
    /// Raised when a machine object tries to access the wrong variant.
    InvalidVariant {
        wrapper: String,
        inner: String,
    },
    /// Raised when wrapper tries to extract an Option<T> from the wrapped type.
    AttributeError {
        obj: String,
        attr: String,
    },
    InvalidInterpolation1dType,
    InvalidInterpolation2dType,
    InvalidPhaseMethod,
    InvalidSteppingMethod,
    InvalidIntersection,
    MachineError(String),
    EvalError(String),
}

impl std::fmt::Display for DexterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PyErr(err) => write!(f, "[D] PyO3 PyErr: '{err}'"),
            Self::CastError(err) => write!(f, "[D] PyO3 CastError: '{err}'"),
            Self::InvalidVariant { wrapper, inner } => write!(
                f,
                "[D] InvalidVariant: '{wrapper}' tried to access non-existent '{inner}' inner type"
            ),
            Self::AttributeError { obj, attr } => write!(
                f,
                "[D] AttributeError: '{obj}' object has no attribute '{attr}'"
            ),
            Self::InvalidInterpolation1dType => write!(
                f,
                concat!(
                    "[D] Supported 1D interpolation types are ",
                    "'Linear', 'Cubic', 'CubicPeriodic', 'Akima', 'AkimaPeriodic' and 'Steffen'",
                )
            ),
            Self::InvalidInterpolation2dType => write!(
                f,
                concat!(
                    "[D] Supported 2D interpolation types are ",
                    "'Bilinear' and 'Bicubic'",
                )
            ),
            Self::InvalidPhaseMethod => write!(
                f,
                concat!(
                    "[D] Supported phase methods are ",
                    "'Zero', 'Average', 'Interpolation' and '('Custom', <value>)'",
                )
            ),
            Self::InvalidSteppingMethod => write!(
                f,
                concat!(
                    "[D] Supported phase methods are ",
                    "'EnergyAdaptiveStep', 'ErrorAdaptiveStep' and '('FixedStep', <value>)'",
                )
            ),
            Self::InvalidIntersection => write!(
                f,
                concat!(
                    "[D] Supported intersection options are ",
                    "'ConstTheta' and 'ConstZeta'",
                )
            ),
            Self::MachineError(err) => write!(f, "[D] MachineError: '{err}'"),
            Self::EvalError(err) => write!(f, "[D] EvalError: '{err}'"),
        }
    }
}

impl From<DexterError> for PyErr {
    fn from(err: DexterError) -> Self {
        PyException::new_err(err.to_string())
    }
}

impl From<PyErr> for DexterError {
    fn from(err: PyErr) -> Self {
        DexterError::PyErr(err.to_string())
    }
}

impl<'a, 'py> From<CastError<'a, 'py>> for DexterError {
    fn from(err: CastError) -> Self {
        DexterError::PyErr(err.to_string())
    }
}

impl From<MachineError> for DexterError {
    fn from(err: MachineError) -> Self {
        DexterError::MachineError(err.to_string())
    }
}

impl From<EvalError> for DexterError {
    fn from(err: EvalError) -> Self {
        DexterError::EvalError(err.to_string())
    }
}
