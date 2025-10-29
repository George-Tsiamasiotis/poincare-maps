mod bfield;
mod current;
mod qfactor;

pub use bfield::PyBfield;
pub use current::PyCurrent;
pub use qfactor::PyQfactor;

/// Generates getter pymethods that return a 1D numpy array.
#[macro_export]
macro_rules! to_numpy1D_impl {
    ($py_object:ident, $rust_object:ident, $name:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $name<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
                self.$rust_object.$name().into_pyarray(py)
            }
        }
    };
}

/// Generates getter pymethods that return a 2D numpy array.
#[macro_export]
macro_rules! to_numpy2D_impl {
    ($py_object:ident, $rust_object:ident, $name:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $name<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
                self.$rust_object.$name().into_pyarray(py)
            }
        }
    };
}

/// Generates a 1D eval method from the wrapped Rust object.
#[macro_export]
macro_rules! eval1D_impl {
    ($py_object:ident, $rust_object:ident, $name:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $name(&mut self, psip: f64) -> Result<f64, PyEqError> {
                Ok(self.$rust_object.$name(psip, &mut self.acc)?)
            }
        }
    };
}

/// Generates a 2D eval method from the wrapped Rust object.
#[macro_export]
macro_rules! eval2D_impl {
    ($py_object:ident, $rust_object:ident, $name:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $name(&mut self, psip: f64, theta: f64) -> Result<f64, PyEqError> {
                Ok(self.$rust_object.$name(
                    psip,
                    theta,
                    &mut self.xacc,
                    &mut self.yacc,
                    &mut self.cache,
                )?)
            }
        }
    };
}

/// Generates a `__repr__` method, corresponding to the wrapped Rust object's `Debug`
/// representation.
#[macro_export]
macro_rules! repr_impl {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn __repr__(&self) -> String {
                format!("{:#?}", self)
            }
        }
    };
}
