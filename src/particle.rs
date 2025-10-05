use crate::point::PointTuple;
use crate::{Bfield, Current, Qfactor};
use crate::{InitialConditions, Point, State};
use crate::{Result, Rk45State};

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pyclass]
#[allow(dead_code)]
pub struct Particle {
    initial: InitialConditions,
    state: State,
    rk45state: Rk45State,
    points: Vec<Point>,
}

#[pymethods]
impl Particle {
    #[new]
    pub fn new(initial: InitialConditions) -> Self {
        Self {
            initial: initial.to_owned(),
            state: State::new(&initial),
            rk45state: Rk45State::default(),
            points: Vec::<Point>::with_capacity(500), // Avoid some relocations?
        }
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

    pub fn get_points(&self) -> Vec<PointTuple> {
        let mut points: Vec<PointTuple> = Vec::with_capacity(self.points.len());
        self.points.iter().for_each(|p| points.push(p.to_tuple()));
        points
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
            self.rk45state.start(&self.state);
            self.rk45state.calculate_k1();
            self.rk45state
                .calculate_state_k2(h, qfactor, bfield, current)?;
            self.rk45state
                .calculate_state_k3(h, qfactor, bfield, current)?;
            self.rk45state
                .calculate_state_k4(h, qfactor, bfield, current)?;
            self.rk45state.calculate_add_terms();
            self.rk45state.calculate_next_state(h);
            self.state = self.rk45state.next.clone();

            self.state.evaluate(qfactor, current, bfield)?;

            self.points.push(Point {
                t: self.state.t,
                theta: self.state.theta,
                psip: self.state.psip,
                rho: self.state.rho,
                zeta: self.state.zeta,
            });
        }
        Ok(())
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
        let initial = InitialConditions {
            t0: 0.0,
            theta0: 0.0,
            psip0: 0.01,
            rho0: 0.00,
            zeta0: 0.1,
            mu: 1e-6,
            pzeta: -0.01,
        };
        let mut particle = Particle::new(initial);
        particle.run(&qfactor, &bfield, &current, 1e-3, 10).unwrap();
    }
}
