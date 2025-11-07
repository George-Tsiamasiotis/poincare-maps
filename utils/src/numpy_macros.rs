//! Wrappers that expose the inner types' [`ArrayD`] methods.
//!
//! The inner type's getter methods must have one of the signatures:
//!     1. `fn(&self) -> Array1<f64>` -> py_get_numpy1D
//!     2. `fn(&self) -> Array2<f64>` -> py_get_numpy2D
//!     1. `fn(&self) -> Result<Array1<f64>, _>` -> py_get_numpy1D_fallible
//!     1. `fn(&self) -> Result<Array2<f64>, _>` -> py_get_numpy2D_fallible
//!
//! It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
//! allow macros to be used inside it.

/// Generates a getter method that returns a 1D numpy array.
#[macro_export]
macro_rules! py_get_numpy1D {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
                self.0.$getter().into_pyarray(py)
            }
        }
    };
}

/// Generates a getter method that returns a 2D numpy array.
#[macro_export]
macro_rules! py_get_numpy2D {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
                self.0.$getter().into_pyarray(py)
            }
        }
    };
}

/// Generates a getter method that returns a 1D numpy array, but the inner type's getter might
/// fail.
#[macro_export]
macro_rules! py_get_numpy1D_fallible {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter<'py>(
                &self,
                py: Python<'py>,
            ) -> Result<Bound<'py, PyArray1<f64>>, PyEqError> {
                Ok(self.0.$getter()?.into_pyarray(py))
            }
        }
    };
}

/// Generates a getter method that returns a 2D numpy array, but the inner type's getter might
/// fail.
#[macro_export]
macro_rules! py_get_numpy2D_fallible {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter<'py>(
                &self,
                py: Python<'py>,
            ) -> Result<Bound<'py, PyArray2<f64>>, PyEqError> {
                Ok(self.0.$getter()?.into_pyarray(py))
            }
        }
    };
}
