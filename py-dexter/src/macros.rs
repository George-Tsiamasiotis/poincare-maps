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
