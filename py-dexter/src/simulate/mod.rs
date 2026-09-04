mod initial;
mod particle;

pub use initial::*;
pub use particle::*;

use dexter::dexter_simulate::{IntersectParams, Intersection, SolverParams};
use pyo3::prelude::*;

use crate::{DexterError, resolve_stepping_method};

// ===============================================================================================

/// This type is only to be used internally when calling integration routines.
#[pyclass(name = "_PySolverParams", frozen, immutable_type)]
pub struct PySolverParams(SolverParams);

#[pymethods]
impl PySolverParams {
    #[new]
    pub fn new<'py>(
        stepping_method: Option<Bound<'py, PyAny>>,
        max_steps: Option<usize>,
        first_step: Option<f64>,
        safety_factor: Option<f64>,
        energy_rel_tol: Option<f64>,
        energy_abs_tol: Option<f64>,
        error_rel_tol: Option<f64>,
        error_abs_tol: Option<f64>,
    ) -> PyResult<Self> {
        let mut solver_params = SolverParams::default();
        if let Some(method) = stepping_method {
            solver_params.method = resolve_stepping_method(method)?;
        }
        max_steps.inspect(|v| solver_params.max_steps = *v);
        first_step.inspect(|v| solver_params.first_step = *v);
        safety_factor.inspect(|v| solver_params.safety_factor = *v);
        energy_rel_tol.inspect(|v| solver_params.energy_rel_tol = *v);
        energy_abs_tol.inspect(|v| solver_params.energy_abs_tol = *v);
        error_rel_tol.inspect(|v| solver_params.error_rel_tol = *v);
        error_abs_tol.inspect(|v| solver_params.error_abs_tol = *v);

        Ok(Self(solver_params))
    }
}

// ===============================================================================================

#[pyclass(name = "_PyIntersectParams", frozen, immutable_type)]
pub struct PyIntersectParams(pub(crate) IntersectParams);

#[pymethods]
impl PyIntersectParams {
    #[new]
    pub fn new<'py>(intersection: String, angle: f64, turns: usize) -> PyResult<Self> {
        let intersection = match intersection.to_lowercase().as_str() {
            "consttheta" => Intersection::ConstTheta,
            "constzeta" => Intersection::ConstZeta,
            _ => return Err(PyErr::from(DexterError::InvalidIntersection)),
        };
        let intersect_params = IntersectParams::new(intersection, angle, turns);

        Ok(Self(intersect_params))
    }
}
