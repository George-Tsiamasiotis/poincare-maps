use equilibrium::*;

use pyo3::prelude::*;

#[pyclass]
#[pyo3(name = "Qfactor")]
pub struct PyQfactor {
    pub qfactor: Qfactor,
}
