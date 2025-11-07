use std::fmt::Debug;
use std::time::Duration;

use utils::array1D_getter_impl;

use crate::State;
use crate::{Distance, Flux, Radians, Time};

use ndarray::Array1;

/// Time series for a Particle's orbit.
#[derive(Clone)]
pub struct Evolution {
    pub time: Vec<Time>,
    pub theta: Vec<Radians>,
    pub psip: Vec<Flux>,
    pub rho: Vec<Distance>,
    pub zeta: Vec<Radians>,
    pub psi: Vec<Flux>,
    pub ptheta: Vec<f64>,
    pub pzeta: Vec<f64>,
    pub duration: Duration,
    pub steps: usize,
}

impl Evolution {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            time: Vec::with_capacity(capacity),
            theta: Vec::with_capacity(capacity),
            psip: Vec::with_capacity(capacity),
            rho: Vec::with_capacity(capacity),
            zeta: Vec::with_capacity(capacity),
            psi: Vec::with_capacity(capacity),
            ptheta: Vec::with_capacity(capacity),
            pzeta: Vec::with_capacity(capacity),
            duration: Duration::default(),
            steps: 0,
        }
    }

    pub fn steps_taken(&self) -> usize {
        self.steps
    }

    pub fn steps_stored(&self) -> usize {
        self.time.len()
    }

    pub fn push_state(&mut self, state: &State) {
        self.time.push(state.time);
        self.theta.push(state.theta);
        self.psip.push(state.psip);
        self.rho.push(state.rho);
        self.zeta.push(state.zeta);
        self.psi.push(state.psi);
        self.ptheta.push(state.ptheta);
        self.pzeta.push(state.pzeta);
    }

    pub fn finish(&mut self) {
        self.time.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.psi.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
    }

    array1D_getter_impl!(time, time, Time);
    array1D_getter_impl!(theta, theta, Radians);
    array1D_getter_impl!(psip, psip, Flux);
    array1D_getter_impl!(rho, rho, Distance);
    array1D_getter_impl!(zeta, zeta, Radians);
    array1D_getter_impl!(psi, psi, Flux);
    array1D_getter_impl!(ptheta, ptheta, f64);
    array1D_getter_impl!(pzeta, pzeta, f64);
}

impl Debug for Evolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Evolution")
            .field(
                "time",
                &format!(
                    "[{:.5}, {:.5}]",
                    self.time.first().unwrap_or(&Time::NAN),
                    self.time.last().unwrap_or(&Time::NAN),
                ),
            )
            .field("duration", &self.duration)
            .field("steps taken", &self.steps_taken())
            .field("steps stored", &self.steps_stored())
            .finish()
    }
}
