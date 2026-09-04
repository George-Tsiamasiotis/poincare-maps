//! Defines the `Particle` wrapper.

use crate::*;
use dexter::dexter_machine::*;
use dexter::dexter_simulate::*;

use pyo3::prelude::*;

#[pyclass(name = "_PyParticle")]
pub struct PyParticle(pub(crate) Particle);

// ===============================================================================================

#[pymethods] // Routines
impl PyParticle {
    #[new]
    pub fn new(initial: &PyInitialConditions) -> Self {
        Self(Particle::new(&initial.0))
    }

    pub fn integrate(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        teval: (f64, f64),
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0.integrate(machine, teval, &solver_params.0);
    }

    pub fn intersect(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        intersect_params: &PyIntersectParams,
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0
            .intersect(machine, &intersect_params.0, &solver_params.0);
    }

    pub fn close(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        periods: usize,
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0.close(machine, periods, &solver_params.0);
    }

    pub fn classify(&mut self, qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        self.0.classify(machine);
    }
}

#[pymethods] // Getters
impl PyParticle {
    #[getter]
    pub fn steps_taken(&self) -> usize {
        self.0.steps_taken()
    }
}

// ===============================================================================================

wrapper_debug_export!(PyParticle);
impl_py_repr!(PyParticle, pretty);
