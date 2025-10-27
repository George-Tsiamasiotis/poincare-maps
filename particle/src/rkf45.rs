use self::tableau::*;
use crate::*;
use equilibrium::{Bfield, Current, Perturbation, Qfactor};

const SAFETY_FACTOR: f64 = 0.9;
const REL_TOL: f64 = 1e-11;

/// Runge-kutta-Fehlberg method coefficients (Wikipedia)
mod tableau {

    // Fifth order
    pub(super) const B1: f64 = 25.0 / 216.0;
    pub(super) const B2: f64 = 0.0;
    pub(super) const B3: f64 = 1408.0 / 2565.0;
    pub(super) const B4: f64 = 2197.0 / 4104.0;
    pub(super) const B5: f64 = -1.0 / 5.0;
    pub(super) const B6: f64 = 0.0;

    // Embedded
    pub(super) const B1E: f64 = 16.0 / 135.0;
    pub(super) const B2E: f64 = 0.0;
    pub(super) const B3E: f64 = 6656.0 / 12825.0;
    pub(super) const B4E: f64 = 28561.0 / 56430.0;
    pub(super) const B5E: f64 = -9.0 / 50.0;
    pub(super) const B6E: f64 = 2.0 / 55.0;

    // pub(super) const C1: f64 = 0.0; Always equals to 0.
    pub(super) const C2: f64 = 1.0 / 4.0;
    pub(super) const C3: f64 = 3.0 / 8.0;
    pub(super) const C4: f64 = 12.0 / 13.0;
    pub(super) const C5: f64 = 1.0;
    pub(super) const C6: f64 = 1.0 / 2.0;

    pub(super) const A21: f64 = 1.0 / 4.0;

    pub(super) const A31: f64 = 3.0 / 32.0;
    pub(super) const A32: f64 = 9.0 / 32.0;

    pub(super) const A41: f64 = 1932.0 / 2197.0;
    pub(super) const A42: f64 = -7200.0 / 2197.0;
    pub(super) const A43: f64 = 7296.0 / 2197.0;

    pub(super) const A51: f64 = 439.0 / 216.0;
    pub(super) const A52: f64 = -8.0;
    pub(super) const A53: f64 = 3680.0 / 513.0;
    pub(super) const A54: f64 = -845.0 / 4104.0;

    pub(super) const A61: f64 = -8.0 / 27.0;
    pub(super) const A62: f64 = 2.0;
    pub(super) const A63: f64 = -3544.0 / 2565.0;
    pub(super) const A64: f64 = 1859.0 / 4104.0;
    pub(super) const A65: f64 = -11.0 / 40.0;
}

pub(crate) struct Solver {
    pub(crate) k1: [f64; 4],
    pub(crate) k2: [f64; 4],
    pub(crate) k3: [f64; 4],
    pub(crate) k4: [f64; 4],
    pub(crate) k5: [f64; 4],
    pub(crate) k6: [f64; 4],
    pub(crate) weights: [f64; 4],
    pub(crate) errors: [f64; 4],
    pub(crate) state1: State,
    pub(crate) state2: State,
    pub(crate) state3: State,
    pub(crate) state4: State,
    pub(crate) state5: State,
    pub(crate) state6: State,
}

impl Solver {
    pub(crate) fn init(&mut self, state: &State) {
        self.state1 = state.clone();
    }

    pub(crate) fn start(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> Result<()> {
        self.calculate_k1();
        self.calculate_state_k2(h, qfactor, bfield, current, per)?;
        self.calculate_state_k3(h, qfactor, bfield, current, per)?;
        self.calculate_state_k4(h, qfactor, bfield, current, per)?;
        self.calculate_state_k5(h, qfactor, bfield, current, per)?;
        self.calculate_state_k6(h, qfactor, bfield, current, per)?;
        self.calculate_embedded_weights();
        self.calculate_errors();
        Ok(())
    }

    pub(crate) fn calculate_k1(&mut self) {
        self.k1 = [
            self.state1.theta_dot,
            self.state1.psip_dot,
            self.state1.rho_dot,
            self.state1.zeta_dot,
        ];
    }

    pub(crate) fn calculate_state_k2(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> Result<()> {
        let coef = [
            A21 * self.k1[0],
            A21 * self.k1[1],
            A21 * self.k1[2],
            A21 * self.k1[3],
        ];

        // States 2-6 are uninitialized, so we need to copy all the indepenendent variables before
        // evaluating.
        self.state2.mu = self.state1.mu;
        self.state2.time = self.state1.time + C2 * h;
        self.state2.theta = self.state1.theta + coef[0] * h;
        self.state2.psip = self.state1.psip + coef[1] * h;
        self.state2.rho = self.state1.rho + coef[2] * h;
        self.state2.zeta = self.state1.zeta + coef[3] * h;
        self.state2.xacc = self.state1.xacc;
        self.state2.yacc = self.state1.yacc;
        self.state2.hcache = self.state1.hcache.clone();
        self.state2.evaluate(qfactor, current, bfield, per)?;
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
        per: &Perturbation,
    ) -> Result<()> {
        let coef = [
            A31 * self.k1[0] + A32 * self.k2[0],
            A31 * self.k1[1] + A32 * self.k2[1],
            A31 * self.k1[2] + A32 * self.k2[2],
            A31 * self.k1[3] + A32 * self.k2[3],
        ];

        self.state3.mu = self.state1.mu;
        self.state3.time = self.state1.time + C3 * h;
        self.state3.theta = self.state1.theta + coef[0] * h;
        self.state3.psip = self.state1.psip + coef[1] * h;
        self.state3.rho = self.state1.rho + coef[2] * h;
        self.state3.zeta = self.state1.zeta + coef[3] * h;
        self.state3.xacc = self.state2.xacc;
        self.state3.yacc = self.state2.yacc;
        self.state3.hcache = self.state2.hcache.clone();
        self.state3.evaluate(qfactor, current, bfield, per)?;
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
        per: &Perturbation,
    ) -> Result<()> {
        let coef = [
            A41 * self.k1[0] + A42 * self.k2[0] + A43 * self.k3[0],
            A41 * self.k1[1] + A42 * self.k2[1] + A43 * self.k3[1],
            A41 * self.k1[2] + A42 * self.k2[2] + A43 * self.k3[2],
            A41 * self.k1[3] + A42 * self.k2[3] + A43 * self.k3[3],
        ];

        self.state4.mu = self.state1.mu;
        self.state4.time = self.state1.time + C4 * h;
        self.state4.theta = self.state1.theta + coef[0] * h;
        self.state4.psip = self.state1.psip + coef[1] * h;
        self.state4.rho = self.state1.rho + coef[2] * h;
        self.state4.zeta = self.state1.zeta + coef[3] * h;
        self.state4.xacc = self.state3.xacc;
        self.state4.yacc = self.state3.yacc;
        self.state4.hcache = self.state3.hcache.clone();
        self.state4.evaluate(qfactor, current, bfield, per)?;
        self.k4 = [
            self.state4.theta_dot,
            self.state4.psip_dot,
            self.state4.rho_dot,
            self.state4.zeta_dot,
        ];
        Ok(())
    }

    pub(crate) fn calculate_state_k5(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> Result<()> {
        let coef = [
            A51 * self.k1[0] + A52 * self.k2[0] + A53 * self.k3[0] + A54 * self.k4[0],
            A51 * self.k1[1] + A52 * self.k2[1] + A53 * self.k3[1] + A54 * self.k4[1],
            A51 * self.k1[2] + A52 * self.k2[2] + A53 * self.k3[2] + A54 * self.k4[2],
            A51 * self.k1[3] + A52 * self.k2[3] + A53 * self.k3[3] + A54 * self.k4[3],
        ];

        self.state5.mu = self.state1.mu;
        self.state5.time = self.state1.time + C5 * h;
        self.state5.theta = self.state1.theta + coef[0] * h;
        self.state5.psip = self.state1.psip + coef[1] * h;
        self.state5.rho = self.state1.rho + coef[2] * h;
        self.state5.zeta = self.state1.zeta + coef[3] * h;
        self.state5.xacc = self.state4.xacc;
        self.state5.yacc = self.state4.yacc;
        self.state5.hcache = self.state4.hcache.clone();
        self.state5.evaluate(qfactor, current, bfield, per)?;
        self.k5 = [
            self.state5.theta_dot,
            self.state5.psip_dot,
            self.state5.rho_dot,
            self.state5.zeta_dot,
        ];
        Ok(())
    }

    pub(crate) fn calculate_state_k6(
        &mut self,
        h: f64,
        qfactor: &Qfactor,
        bfield: &Bfield,
        current: &Current,
        per: &Perturbation,
    ) -> Result<()> {
        #[rustfmt::skip]
            let coef = [
                A61 * self.k1[0] + A62 * self.k2[0] + A63 * self.k3[0] + A64 * self.k4[0] + A65 * self.k5[0],
                A61 * self.k1[1] + A62 * self.k2[1] + A63 * self.k3[1] + A64 * self.k4[1] + A65 * self.k5[1],
                A61 * self.k1[2] + A62 * self.k2[2] + A63 * self.k3[2] + A64 * self.k4[2] + A65 * self.k5[2],
                A61 * self.k1[3] + A62 * self.k2[3] + A63 * self.k3[3] + A64 * self.k4[3] + A65 * self.k5[3],
            ];

        self.state6.mu = self.state1.mu;
        self.state6.time = self.state1.time + C6 * h;
        self.state6.theta = self.state1.theta + coef[0] * h;
        self.state6.psip = self.state1.psip + coef[1] * h;
        self.state6.rho = self.state1.rho + coef[2] * h;
        self.state6.zeta = self.state1.zeta + coef[3] * h;
        self.state6.xacc = self.state5.xacc;
        self.state6.yacc = self.state5.yacc;
        self.state6.hcache = self.state5.hcache.clone();
        self.state6.evaluate(qfactor, current, bfield, per)?;
        self.k6 = [
            self.state6.theta_dot,
            self.state6.psip_dot,
            self.state6.rho_dot,
            self.state6.zeta_dot,
        ];
        Ok(())
    }

    pub(crate) fn calculate_embedded_weights(&mut self) {
        for i in 0..self.weights.len() {
            self.weights[i] = B1 * self.k1[i]
                + B2E * self.k2[i]
                + B3E * self.k3[i]
                + B4E * self.k4[i]
                + B5E * self.k5[i]
                + B6E * self.k6[i];
        }
    }
    pub(crate) fn calculate_errors(&mut self) {
        for i in 0..self.errors.len() {
            self.errors[i] = (B1 - B1E) * self.k1[i]
                + (B2 - B2E) * self.k2[i]
                + (B3 - B3E) * self.k3[i]
                + (B4 - B4E) * self.k4[i]
                + (B5 - B5E) * self.k5[i]
                + (B6 - B6E) * self.k6[i];
        }
    }

    /// Adjust the error by calculating the relative difference in the energy at every step.
    #[cfg(not(feature = "error-adaptive-step"))]
    pub(crate) fn calculate_optimal_step(&mut self, h: f64) -> f64 {
        let initial_energy = self.state1.energy();
        let final_energy = self.state6.energy();
        // // When the energy happens to be smaller than REL_TOL, the optimal step keeps getting
        // // smaller due to the `REL_TOL/energy_diff` factor, so we need to bound it
        let energy_diff = ((initial_energy - final_energy) / initial_energy)
            .abs()
            .max(REL_TOL / 2.0);
        let exp = if energy_diff >= REL_TOL { 0.2 } else { 0.25 };
        SAFETY_FACTOR * h * (REL_TOL / energy_diff).powf(exp)
    }

    /// Source:
    /// https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ss2017/numerische_Methoden_fuer_komplexe_Systeme_II/rkm-1.pdf
    #[cfg(feature = "error-adaptive-step")]
    pub(crate) fn calculate_optimal_step(&mut self, h: f64) -> f64 {
        // Using the max error vs each variable's error is equivalent.

        let mut max_error = (*self
            .errors
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap())
        .abs();
        // When all errors happen to be smaller than REL_TOL, the optimal step keeps getting
        // smaller due to the `REL_TOL/max_error` factor, so we need to bound it
        max_error = max_error.max(REL_TOL / 2.0);

        // 0.2 = 1/*(p+1), where p the order
        let exp = if max_error >= REL_TOL { 0.2 } else { 0.25 };
        SAFETY_FACTOR * h * (REL_TOL / max_error).powf(exp)
    }

    pub(crate) fn next_state(&mut self, h: f64) -> State {
        {
            State {
                time: self.state1.time + h,
                theta: self.state1.theta + h * self.weights[0],
                psip: self.state1.psip + h * self.weights[1],
                rho: self.state1.rho + h * self.weights[2],
                zeta: self.state1.zeta + h * self.weights[3],
                mu: self.state1.mu,
                xacc: self.state6.xacc,
                yacc: self.state6.yacc,
                hcache: self.state6.hcache.clone(),
                ..Default::default()
            }
        }
    }
}

impl Default for Solver {
    fn default() -> Self {
        Self {
            k1: [f64::NAN; 4],
            k2: [f64::NAN; 4],
            k3: [f64::NAN; 4],
            k4: [f64::NAN; 4],
            k5: [f64::NAN; 4],
            k6: [f64::NAN; 4],
            weights: [f64::NAN; 4],
            errors: [f64::NAN; 4],
            state1: State::default(),
            state2: State::default(),
            state3: State::default(),
            state4: State::default(),
            state5: State::default(),
            state6: State::default(),
        }
    }
}
