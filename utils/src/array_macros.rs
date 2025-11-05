/// Generates getters that return \[T\] fields to Array1<T>
///
/// This is needed for implementing python getter wrappers.
#[macro_export]
macro_rules! array1D_getter_impl {
    ($object:ident, $fun_name:ident, $vec_name:ident $(.$field:ident),*) => {
        impl $object {
            pub fn $fun_name(&self) -> Array1<f64> {
                Array1::from(self.$vec_name$(.$field),*.clone())
            }
        }
    };
}

/// Generates getters that return 2D arrays.
///
/// This is needed for implementing python getter wrappers.
#[macro_export]
macro_rules! array2D_getter_impl {
    ($object:ident, $fun_name:ident, $vec_name:ident $(.$field:ident),*) => {
        impl $object {
            pub fn $fun_name(&self) -> Array2<f64> {
                Array2::from(self.$vec_name$(.$field),*.clone())
            }
        }
    };
}
