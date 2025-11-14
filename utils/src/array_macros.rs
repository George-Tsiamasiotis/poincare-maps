/// Generates getters that return `[T]` fields to `Array1<T>`.
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
#[macro_export]
macro_rules! array2D_getter_impl {
    // field
    ($fun_name:ident, $($field:ident).+, $return_type:ty) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($($field).+)]
        #[doc = "` field as an [`Array2<"]
        #[doc = stringify!($return_type>)]
        #[doc = "`]" ]
        pub fn $fun_name(&self) -> Array2<$return_type> {
            Array2::from(self.$($field).+.clone())
        }
    };

    // method
    ($fun_name:ident, $($field:ident).+, (), $return_type:ty) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($($field).+)]
        #[doc = "` field as an [`Array2<"]
        #[doc = stringify!($return_type>)]
        #[doc = "`]" ]
        pub fn $fun_name(&self) -> Array2<$return_type> {
            Array2::from(self.$($field).+().clone())
        }
    }
}
