//! Miscellaneous getter methods
//!
//! It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
//! allow macros to be used inside it.

/// Generates a impl Debug block, using the innet type's Debug representation.
#[macro_export]
macro_rules! py_debug_impl {
    ($py_object:ident) => {
        impl std::fmt::Debug for $py_object {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                self.0.fmt(f)
            }
        }
    };
}

/// Generates a `__repr__` method, corresponding to the inner type's `Debug` representation.
#[macro_export]
macro_rules! py_repr_impl {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn __repr__(&self) -> String {
                format!("{:#?}", self.0)
            }
        }
    };
}

/// Generates a getter method for a float.
///
/// The float should be returned by the wrapped type as `py_object.0.<field_name>()`.
#[macro_export]
macro_rules! py_get_float {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> f64 {
                self.0.$getter()
            }
        }
    };
}

/// Generates a getter method for the NetCDF's path field.
#[macro_export]
macro_rules! py_get_path {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn path(&self) -> String {
                String::from(safe_unwrap!(
                    "file already opened",
                    self.0.path.clone().into_os_string().into_string()
                ))
            }
        }
    };
}

/// Generates a getter method for the object's interpolation type.
#[macro_export]
macro_rules! py_get_typ {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn typ(&self) -> String {
                String::from(self.0.typ.clone())
            }
        }
    };
}
