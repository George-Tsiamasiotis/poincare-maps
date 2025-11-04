use std::fmt::Debug;
use std::time::Duration;

use ndarray::Array1;
use utils::array1D_getter_impl;

use crate::Point;

/// Time series for a Particle's orbit.
#[derive(Clone)]
pub struct Evolution {
    pub time: Vec<f64>,
    pub theta: Vec<f64>,
    pub psip: Vec<f64>,
    pub rho: Vec<f64>,
    pub zeta: Vec<f64>,
    pub psi: Vec<f64>,
    pub ptheta: Vec<f64>,
    pub pzeta: Vec<f64>,
    pub duration: Duration,
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
        }
    }

    pub fn push_point(&mut self, point: &Point) {
        self.time.push(point.time);
        self.theta.push(point.theta);
        self.psip.push(point.psip);
        self.rho.push(point.rho);
        self.zeta.push(point.zeta);
        self.psi.push(point.psi);
        self.ptheta.push(point.ptheta);
        self.pzeta.push(point.pzeta);
    }

    pub fn shrink_to_fit(&mut self) {
        self.time.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.psi.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
    }
}

array1D_getter_impl!(Evolution, time, time);
array1D_getter_impl!(Evolution, theta, theta);
array1D_getter_impl!(Evolution, psip, psip);
array1D_getter_impl!(Evolution, rho, rho);
array1D_getter_impl!(Evolution, zeta, zeta);
array1D_getter_impl!(Evolution, psi, psi);
array1D_getter_impl!(Evolution, ptheta, ptheta);
array1D_getter_impl!(Evolution, pzeta, pzeta);

impl Debug for Evolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Evolution")
            .field(
                "time",
                &format!(
                    "[{:.5}, {:.5}]",
                    self.time.first().unwrap_or(&f64::NAN),
                    self.time.last().unwrap_or(&f64::NAN),
                ),
            )
            .field("duration", &self.duration)
            .field("length", &self.time.len())
            .finish()
    }
}
