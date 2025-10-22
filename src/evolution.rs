use core::f64;
use std::time::Duration;

use crate::EVOLUTION_INIT_CAPACITY;
use crate::Point;
use crate::State;
use crate::numpy_getter_1D;

use numpy::{PyArray1, ToPyArray};
use pyo3::prelude::*;

/// Time series of a particle's orbit.
#[derive(Debug, Clone)]
#[pyclass]
pub struct Evolution {
    pub time: Vec<f64>,
    pub theta: Vec<f64>,
    pub psip: Vec<f64>,
    pub rho: Vec<f64>,
    pub zeta: Vec<f64>,
    pub ptheta: Vec<f64>,
    pub pzeta: Vec<f64>,
    pub psi: Vec<f64>,
    pub mu: Vec<f64>, // Constant for now, but might change
    pub duration: Duration,
    pub steps_taken: usize,
    pub steps_stored: usize,
}

numpy_getter_1D!(Evolution, time);
numpy_getter_1D!(Evolution, theta);
numpy_getter_1D!(Evolution, psip);
numpy_getter_1D!(Evolution, rho);
numpy_getter_1D!(Evolution, zeta);
numpy_getter_1D!(Evolution, ptheta);
numpy_getter_1D!(Evolution, psi);
numpy_getter_1D!(Evolution, mu);

impl Evolution {
    pub fn new() -> Self {
        Self::default()
    }

    /// Initializes the time series vectors with a specified capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            time: Vec::with_capacity(capacity),
            theta: Vec::with_capacity(capacity),
            psip: Vec::with_capacity(capacity),
            rho: Vec::with_capacity(capacity),
            zeta: Vec::with_capacity(capacity),
            ptheta: Vec::with_capacity(capacity),
            pzeta: Vec::with_capacity(capacity),
            psi: Vec::with_capacity(capacity),
            mu: Vec::with_capacity(capacity),
            duration: Duration::default(),
            steps_taken: 0,
            steps_stored: 0,
        }
    }

    /// Creates a non-evaluated [`State`] from the most recently pushed values.
    pub fn last_state(&self) -> State {
        State::from_initial(
            self.time.last().copied().unwrap_or(f64::NAN),
            self.theta.last().copied().unwrap_or(f64::NAN),
            self.psip.last().copied().unwrap_or(f64::NAN),
            self.rho.last().copied().unwrap_or(f64::NAN),
            self.zeta.last().copied().unwrap_or(f64::NAN),
            self.mu.last().copied().unwrap_or(f64::NAN),
        )
    }

    /// Adds a [`Point`] to the time series vector.
    pub fn push_point(&mut self, point: Point) {
        self.time.push(point.time);
        self.theta.push(point.theta);
        self.psip.push(point.psip);
        self.rho.push(point.rho);
        self.zeta.push(point.zeta);
        self.ptheta.push(point.ptheta);
        self.pzeta.push(point.pzeta);
        self.psi.push(point.psi);
    }

    pub fn shrink_to_fit(&mut self) {
        self.time.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
        self.psi.shrink_to_fit();
    }
}

impl Default for Evolution {
    fn default() -> Self {
        Self::with_capacity(EVOLUTION_INIT_CAPACITY)
    }
}
