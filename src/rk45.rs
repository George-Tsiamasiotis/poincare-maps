use crate::Bfield;
use crate::Current;
use crate::Qfactor;
use crate::Result;
use crate::State;

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
        self.state2.t = self.state1.t + h / 2.0;
        self.state2.theta = self.state1.theta + h / 2.0 * self.k1[0];
        self.state2.psip = self.state1.psip + h / 2.0 * self.k1[1];
        self.state2.rho = self.state1.rho + h / 2.0 * self.k1[2];
        self.state2.zeta = self.state1.zeta + h / 2.0 * self.k1[3];
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
        self.state3.t = self.state1.t + h / 2.0;
        self.state3.theta = self.state1.theta + h / 2.0 * self.k2[0];
        self.state3.psip = self.state1.psip + h / 2.0 * self.k2[1];
        self.state3.rho = self.state1.rho + h / 2.0 * self.k2[2];
        self.state3.zeta = self.state1.zeta + h / 2.0 * self.k2[3];
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
        self.state4.t = self.state1.t + h;
        self.state4.theta = self.state1.theta + h * self.k3[0];
        self.state4.psip = self.state1.psip + h * self.k3[1];
        self.state4.rho = self.state1.rho + h * self.k3[2];
        self.state4.zeta = self.state1.zeta + h * self.k3[3];
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
            self.weights[i] = (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]) / 6.0;
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
