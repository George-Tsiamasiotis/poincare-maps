use crate::Bfield;
use crate::Current;
use crate::Qfactor;
use crate::Result;
use crate::State;

// Hardcoded. Change this to switch between methods.
use self::classic_tableau::*;

#[derive(Debug)]
pub(crate) struct Rk45State {
    pub(crate) k1: [f64; 4],
    pub(crate) k2: [f64; 4],
    pub(crate) k3: [f64; 4],
    pub(crate) k4: [f64; 4],
    pub(crate) weights: [f64; 4],
    pub(crate) state1: State,
    pub(crate) state2: State,
    pub(crate) state3: State,
    pub(crate) state4: State,
    pub(crate) next: State,
}

impl Rk45State {
    pub(crate) fn init(&mut self, state: &State) {
        // Unnecessary copy makes makes the code a bit easier to follow, will probably be removed.
        // Although maybe its necessary to keep the state of the Accelerator.
        self.state1 = state.to_owned();
        self.state2 = state.to_owned();
        self.state3 = state.to_owned();
        self.state4 = state.to_owned();
        self.next = State::default();
    }

    pub(crate) fn start(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
    ) -> Result<()> {
        self.calculate_k1();
        self.calculate_state_k2(h, qfactor, bfield, current)?;
        self.calculate_state_k3(h, qfactor, bfield, current)?;
        self.calculate_state_k4(h, qfactor, bfield, current)?;
        self.calculate_weights();
        self.calculate_next_state(h);
        Ok(())
    }

    pub(crate) fn calculate_k1(&mut self) {
        self.k1 = [
            self.state1.theta_dot,
            self.state1.psip_dot,
            self.state1.rho_dot,
            self.state1.zeta_dot,
        ]
    }

    pub(crate) fn calculate_state_k2(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
    ) -> Result<()> {
        let coef = [
            A21 * self.k1[0],
            A21 * self.k1[1],
            A21 * self.k1[2],
            A21 * self.k1[3],
        ];

        self.state2.t = self.state1.t + C2 * h;
        self.state2.theta = self.state1.theta + coef[0] * h;
        self.state2.psip = self.state1.psip + coef[1] * h;
        self.state2.rho = self.state1.rho + coef[2] * h;
        self.state2.zeta = self.state1.zeta + coef[3] * h;
        self.state2.evaluate(qfactor, current, bfield)?;
        self.k2 = [
            self.state2.theta_dot,
            self.state2.psip_dot,
            self.state2.rho_dot,
            self.state2.zeta_dot,
        ];
        Ok(())
    }

    pub(crate) fn calculate_state_k3(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
    ) -> Result<()> {
        let coef = [
            A31 * self.k1[0] + A32 * self.k2[0],
            A31 * self.k1[1] + A32 * self.k2[1],
            A31 * self.k1[2] + A32 * self.k2[2],
            A31 * self.k1[3] + A32 * self.k2[3],
        ];
        self.state3.t = self.state1.t + C3 * h;
        self.state3.theta = self.state1.theta + coef[0] * h;
        self.state3.psip = self.state1.psip + coef[1] * h;
        self.state3.rho = self.state1.rho + coef[2] * h;
        self.state3.zeta = self.state1.zeta + coef[3] * h;
        self.state3.evaluate(qfactor, current, bfield)?;
        self.k3 = [
            self.state3.theta_dot,
            self.state3.psip_dot,
            self.state3.rho_dot,
            self.state3.zeta_dot,
        ];
        Ok(())
    }

    pub(crate) fn calculate_state_k4(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
    ) -> Result<()> {
        let coef = [
            A41 * self.k1[0] + A42 * self.k2[0] + A43 * self.k3[0],
            A41 * self.k1[1] + A42 * self.k2[1] + A43 * self.k3[1],
            A41 * self.k1[2] + A42 * self.k2[2] + A43 * self.k3[2],
            A41 * self.k1[3] + A42 * self.k2[3] + A43 * self.k3[3],
        ];
        self.state4.t = self.state1.t + C4 * h;
        self.state4.theta = self.state1.theta + coef[0] * h;
        self.state4.psip = self.state1.psip + coef[1] * h;
        self.state4.rho = self.state1.rho + coef[2] * h;
        self.state4.zeta = self.state1.zeta + coef[3] * h;
        self.state4.evaluate(qfactor, current, bfield)?;
        self.k4 = [
            self.state4.theta_dot,
            self.state4.psip_dot,
            self.state4.rho_dot,
            self.state4.zeta_dot,
        ];
        Ok(())
    }

    pub(crate) fn calculate_weights(&mut self) {
        for i in 0..self.weights.len() {
            self.weights[i] = B1 * self.k1[i] + B2 * self.k2[i] + B3 * self.k3[i] + B4 * self.k4[i];
        }
    }

    pub(crate) fn calculate_next_state(&mut self, h: f64) {
        self.next.t = self.state1.t + h;
        self.next.theta = self.state1.theta + h * self.weights[0];
        self.next.psip = self.state1.psip + h * self.weights[1];
        self.next.rho = self.state1.rho + h * self.weights[2];
        self.next.zeta = self.state1.zeta + h * self.weights[3];
        self.next.mu = self.state1.mu;
    }

    pub(crate) fn next_state(&self) -> State {
        self.next.clone()
    }
}

impl Default for Rk45State {
    fn default() -> Self {
        Self {
            k1: [f64::NAN; 4],
            k2: [f64::NAN; 4],
            k3: [f64::NAN; 4],
            k4: [f64::NAN; 4],
            weights: [f64::NAN; 4],
            state1: State::new_uninit(),
            state2: State::new_uninit(),
            state3: State::new_uninit(),
            state4: State::new_uninit(),
            next: State::new_uninit(),
        }
    }
}

mod classic_tableau {

    pub(super) const B1: f64 = 1.0 / 6.0;
    pub(super) const B2: f64 = 1.0 / 3.0;
    pub(super) const B3: f64 = 1.0 / 3.0;
    pub(super) const B4: f64 = 1.0 / 6.0;

    // pub(super) const C1: f64 = 0.0; Always equals to 0.
    pub(super) const C2: f64 = 1.0 / 2.0;
    pub(super) const C3: f64 = 1.0 / 2.0;
    pub(super) const C4: f64 = 1.0;

    pub(super) const A21: f64 = 1.0 / 2.0;

    pub(super) const A31: f64 = 0.0;
    pub(super) const A32: f64 = 1.0 / 2.0;

    pub(super) const A41: f64 = 0.0;
    pub(super) const A42: f64 = 0.0;
    pub(super) const A43: f64 = 1.0;
}
