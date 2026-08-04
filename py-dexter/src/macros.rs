#[macro_export]
macro_rules! impl_py_repr {
    ($obj: ident, simple) => {
        #[pymethods]
        impl $obj {
            pub fn __repr__(&self) -> String {
                format!("{:?}", self.0)
            }
        }
    };
    ($obj: ident, pretty) => {
        #[pymethods]
        impl $obj {
            pub fn __repr__(&self) -> String {
                format!("{:#?}", self.0)
            }
        }
    };
}

#[macro_export]
macro_rules! wrapper_debug_export {
    ($obj: ident) => {
        impl std::fmt::Debug for $obj {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                self.0.fmt(f)
            }
        }
    };
}
