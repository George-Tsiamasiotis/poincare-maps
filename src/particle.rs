use crate::{Bfield, Current, Qfactor};
use crate::{InitialConditions, State};
use crate::{Result, Rk45State};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Initial capacity of the vectors that store the evolution of the particle.
const VEC_INIT_CAPACITY: usize = 500;

#[pyclass]
pub struct Particle {
    /// The initial (θ, ψ_p, ρ, ζ, μ) of the particle.
    initial: InitialConditions,
    /// The current state of the particle.
    state: State,
    /// The current RK45 state.
    rk45state: Rk45State,
    /// The calculated evaluation times.
    #[pyo3(get)]
    t: Vec<f64>,
    /// The calculated θ values.
    #[pyo3(get)]
    theta: Vec<f64>,
    /// The calculated ψ_p values.
    #[pyo3(get)]
    psip: Vec<f64>,
    /// The calculated ρ values.
    #[pyo3(get)]
    rho: Vec<f64>,
    /// The calculated ζ values.
    #[pyo3(get)]
    zeta: Vec<f64>,
    /// The calculated Pζ values.
    #[pyo3(get)]
    pzeta: Vec<f64>,
    /// The calculated Pθ values.
    #[pyo3(get)]
    ptheta: Vec<f64>,
}

#[pymethods]
impl Particle {
    #[new]
    pub fn new(initial: InitialConditions) -> Self {
        Self {
            initial: initial.to_owned(),
            state: State::new_init(&initial),
            rk45state: Rk45State::default(),
            t: Vec::with_capacity(VEC_INIT_CAPACITY),
            theta: Vec::with_capacity(VEC_INIT_CAPACITY),
            psip: Vec::with_capacity(VEC_INIT_CAPACITY),
            rho: Vec::with_capacity(VEC_INIT_CAPACITY),
            zeta: Vec::with_capacity(VEC_INIT_CAPACITY),
            pzeta: Vec::with_capacity(VEC_INIT_CAPACITY),
            ptheta: Vec::with_capacity(VEC_INIT_CAPACITY),
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }

    #[pyo3(name = "run")]
    pub fn run_py(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        step_size: f64,
        steps: usize,
    ) -> PyResult<()> {
        match self.run(qfactor, bfield, current, step_size, steps) {
            Ok(()) => Ok(()),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }
}

impl Particle {
    pub(crate) fn run(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        step_size: f64,
        steps: usize,
    ) -> Result<()> {
        self.state.evaluate(qfactor, current, bfield)?;

        let h = step_size;
        // TEMP
        for _ in 0..steps {
            self.rk45state = Rk45State::default();
            self.rk45state.init(&self.state);
            self.rk45state.start(h, qfactor, bfield, current)?;
            self.state = self.rk45state.next_state();

            self.state.evaluate(qfactor, current, bfield)?;

            self.t.push(self.state.theta);
            self.theta.push(self.state.theta);
            self.psip.push(self.state.psip);
            self.rho.push(self.state.rho);
            self.zeta.push(self.state.zeta);
            self.pzeta.push(self.state.pzeta);
            self.ptheta.push(self.state.ptheta);
        }

        Ok(())
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let evolution = if self.t.is_empty() {
            "not calculated".to_string()
        } else {
            format!("{} steps", self.t.len())
        };

        f.debug_struct("Particle")
            .field("initial", &self.initial)
            .field("Evolution", &evolution)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_particle() {
        let path = PathBuf::from("./data.nc");
        let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
        let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
        let current = Current::from_dataset(&path, "akima").unwrap();

        // Normal passing particle
        let initial = InitialConditions {
            t0: 0.0,
            theta0: 3.14,
            psip0: 0.05,
            rho0: 0.05,
            zeta0: 0.1,
            mu: 1e-4,
        };
        let mut particle = Particle::new(initial);
        particle
            .run(&qfactor, &bfield, &current, 1e-1, 5000)
            .unwrap();
    }
}
