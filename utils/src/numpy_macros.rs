/// Generates getter pymethods that return a 1D numpy array.
///
/// It is necessary to create a new #[pymethods] impl block every time, since #[pymethods] does not
/// allow macros to be used inside it.
#[macro_export]
macro_rules! to_numpy1D_impl {
    ($py_object:ident, $rust_object:ident, $name:ident $(.$field:ident),*) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $name<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
                self.$rust_object.$name$(.$field),*().into_pyarray(py)
            }
        }
    };
}

/// Generates getter pymethods that return a 2D numpy array.
///
/// It is necessary to create a new #[pymethods] impl block every time, since #[pymethods] does not
/// allow macros to be used inside it.
#[macro_export]
macro_rules! to_numpy2D_impl {
    ($py_object:ident, $rust_object:ident, $name:ident $(.$field:ident),*) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $name<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
                self.$rust_object.$name$(.$field),*().into_pyarray(py)
            }
        }
    };
}
