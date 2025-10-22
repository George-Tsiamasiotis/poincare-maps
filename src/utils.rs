use pyo3::prelude::*;
use pyo3::types::PyString;

/// Defines a getter function exposed to python for a field that can be converted to a 1D Numpy array.
#[macro_export]
macro_rules! numpy_getter_1D {
    ($object: ident, $field: ident) => {
        #[pymethods]
        impl $object {
            #[getter]
            #[coverage(off)]
            pub fn $field<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
                self.$field.to_pyarray(py)
            }
        }
    };
}

/// Defines a getter function exposed to python for a field that can be converted to a 2D Numpy array.
#[macro_export]
macro_rules! numpy_getter_2D {
    ($object: ident, $field: ident) => {
        #[pymethods]
        impl $object {
            #[getter]
            #[coverage(off)]
            pub fn $field<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
                self.$field.to_pyarray(py)
            }
        }
    };
}

#[macro_export]
macro_rules! __repr__ {
    ($object: ident) => {
        #[pymethods]
        impl $object {
            #[coverage(off)]
            pub fn __repr__(&self) -> String {
                format!("{:#?}", &self)
            }
        }
    };
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The particle's integration status
#[derive(Default, Debug, Clone)]
pub enum IntegrationStatus {
    #[default]
    Initialized,
    Integrated,
    Escaped,
    TimedOut,
}

impl<'py> IntoPyObject<'py> for IntegrationStatus {
    type Target = PyString;
    type Output = Bound<'py, Self::Target>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> std::result::Result<Self::Output, Self::Error> {
        use IntegrationStatus::*;
        match self {
            Initialized => Ok("Initialized".into_pyobject(py)?),
            Integrated => Ok("Integrated".into_pyobject(py)?),
            Escaped => Ok("Escaped".into_pyobject(py)?),
            TimedOut => Ok("TimedOut".into_pyobject(py)?),
        }
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Defines the `angle=const` cross section.
#[derive(Debug, Clone)]
pub enum Surface {
    ConstZeta,
    ConstTheta,
}

#[pyclass]
#[derive(Debug, Clone)]
pub struct PoincareParameters {
    pub surface: Surface,
    pub intersection: f64,
    pub turns: usize,
}

#[pymethods]
impl PoincareParameters {
    #[new]
    pub fn new(surface: &str, intersection: f64, turns: usize) -> Self {
        let surface = match surface.to_lowercase().as_str() {
            "zeta" => Surface::ConstZeta,
            "theta" => Surface::ConstTheta,
            _ => panic!("Intersection angle must be either 'theta' or 'zeta'"),
        };
        Self {
            surface,
            intersection,
            turns,
        }
    }
}
