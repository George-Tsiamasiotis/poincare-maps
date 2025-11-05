/// Generates a `__repr__` method, corresponding to the wrapped Rust object's `Debug` representation.
///
/// It is necessary to create a new #[pymethods] impl block every time, since #[pymethods] does not
/// allow macros to be used inside it.
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

#[macro_export]
macro_rules! to_pyfloat_impl {
    ($py_object:ident, $rust_object:ident, $name:ident$(.$field:ident),*) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $name(&self) -> f64 {
                self.$rust_object.$name$(.$field),*.into()
            }
        }
    };
}

/// Generates a getter method for String types.
///
/// Pathbufs need seperate treatment
#[macro_export]
macro_rules! to_pystr_impl {
    ($py_object:ident, $rust_object:ident, path) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn get_path(&self) -> String {
                String::from(
                    self.$rust_object
                        .path
                        .clone()
                        .into_os_string()
                        .into_string()
                        .unwrap_or("".into()),
                )
            }
        }
    };
    ($py_object:ident, $rust_object:ident, $name:ident$(.$field:ident),*) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $name(&self) -> String {
                String::from(self.$rust_object.$name$(.$field),*.clone())
            }
        }
    };
}
