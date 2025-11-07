use equilibrium::EqError;
use particle::ParticleError;
use poincare::PoincareError;
use pyo3::exceptions::PyException;
use pyo3::PyErr;

/// Creates newtype wrappers around the foreign error types, to allow conversion to [`PyErr`] and
/// use of the `?` operator.
///
/// Source:
/// [`https://pyo3.rs/main/function/error-handling#foreign-rust-error-types`]
macro_rules! to_pyerr_impl {
    ($error_type:ident, $py_error_type: ident) => {
        #[derive(Debug)]
        pub struct $py_error_type($error_type);

        impl From<$py_error_type> for PyErr {
            fn from(error: $py_error_type) -> Self {
                PyException::new_err(error.0.to_string())
            }
        }

        impl From<$error_type> for $py_error_type {
            fn from(other: $error_type) -> Self {
                Self(other)
            }
        }
    };
}

to_pyerr_impl!(EqError, PyEqError);
to_pyerr_impl!(ParticleError, PyParticleError);
to_pyerr_impl!(PoincareError, PyPoincareError);
