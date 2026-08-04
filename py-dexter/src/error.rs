use crate::*;
use pyo3::{exceptions::PyException, prelude::*};

#[derive(Debug)]
pub enum DexterError {
    PyErr(String),
    /// Raised when an equilibrium object tries to access the wrong variant.
    InvalidVariant {
        wrapper: String,
        inner: String,
    },
    InvalidInterpolationType(String),
    EqError(String),
    EvalError(String),
}

impl std::fmt::Display for DexterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidVariant { wrapper, inner } => write!(
                f,
                "InvalidVariant: '{wrapper}' tried to access non-existent '{inner}' inner type"
            ),
            _ => write!(f, "{:?}", self),
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

impl From<EqError> for DexterError {
    fn from(err: EqError) -> Self {
        DexterError::EqError(err.to_string())
    }
}

impl From<EvalError> for DexterError {
    fn from(err: EvalError) -> Self {
        DexterError::EvalError(err.to_string())
    }
}
