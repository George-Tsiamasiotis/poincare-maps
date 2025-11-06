/// Generates a 1D eval method from the wrapped Rust object.
///
/// It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
/// allow macros to be used inside it.
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
///
/// It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
/// allow macros to be used inside it.
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

/// Generates an eval method from the wrapped Rust Harmonic object.
///
/// It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
/// allow macros to be used inside it.
#[macro_export]
macro_rules! eval_harmonic_impl {
    ($py_object:ident, $rust_object:ident, $name:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $name(&mut self, psip: f64, theta: f64, zeta: f64) -> Result<f64, PyEqError> {
                Ok(self
                    .$rust_object
                    .$name(psip, theta, zeta, &mut self.cache, &mut self.acc)?)
            }
        }
    };
}
