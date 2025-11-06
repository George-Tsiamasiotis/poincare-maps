/// Generates getters that return `[T]` fields to `Array1<T>`.
///
/// This is needed for implementing python getter wrappers.
#[macro_export]
macro_rules! array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+, $return_type:ty) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($($field).+)]
        #[doc = "` field as an [`Array1<"]
        #[doc = stringify!($return_type>)]
        #[doc = "`]" ]
        pub fn $fun_name(&self) -> Array1<$return_type> {
            Array1::from(self.$($field).+.clone())
        }
    }
}

/// Generates getters that return 2D arrays.
///
/// This is needed for implementing python getter wrappers.
#[macro_export]
macro_rules! array2D_getter_impl {
    ($object:ident, $fun_name:ident, $first_field:ident $(.$field:ident),*) => {
        impl $object {
            #[doc = "returns the `" ]
            #[doc = stringify!($first_field$(.$field),*)]
            #[doc = "` field as a 2D array." ]
            pub fn $fun_name(&self) -> Array2<f64> {
                Array2::from(self.$first_field$(.$field),*.clone())
            }
        }
    };
}
