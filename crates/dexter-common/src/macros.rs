//! Utility macros.

/// Generates a getter method that returns a `Vec<T>` field as an `Array1<T>`.
#[macro_export]
macro_rules! array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values as a 1D array." ]
        pub fn $fun_name(&self) -> Array1<f64> {
            Array1::from(self.$($field).+.clone())
        }
    }
}

/// Generates a getter method for an `Array1<T>` field.
#[macro_export]
macro_rules! export_array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values as a 1D array." ]
        pub fn $fun_name(&self) -> Array1<f64> {
            self.$($field).+.$fun_name()
        }
    }
}
