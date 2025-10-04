use tokamak_equilibria::{Bfield, Current, Qfactor};

use crate::InitialConditions;
use crate::Point;
use crate::Result;
use crate::Rk45State;
use crate::State;

#[allow(dead_code)]
pub struct System<'a> {
    pub(crate) qfactor: &'a Qfactor,
    pub(crate) bfield: &'a Bfield,
    pub(crate) current: &'a Current,
    pub(crate) initial_state: State,
    pub(crate) rk45_state: Rk45State,
    pub(crate) state: State,
    pub(crate) state2: State,
    pub(crate) state3: State,
    pub(crate) state4: State,
    pub(crate) next: State,
    pub(crate) points: Vec<Point>,
}

/// Evaluation methods.
impl<'a> System<'a> {
    pub fn init(
        qfactor: &'a Qfactor,
        current: &'a Current,
        bfield: &'a Bfield,
        initial: InitialConditions,
    ) -> Result<Self> {
        let mut initial_state = State::init(&initial)?;
        initial_state.evaluate(qfactor, current, bfield)?;

        let points: Vec<Point> = Vec::with_capacity(100);

        Ok(Self {
            qfactor,
            bfield,
            current,
            initial_state: initial_state.clone(),
            rk45_state: Rk45State::default(),
            state: initial_state.clone(),
            state2: State::init(&initial)?,
            state3: State::init(&initial)?,
            state4: State::init(&initial)?,
            next: State::init(&initial)?,
            points,
        })
    }

    #[allow(dead_code)]
    fn step(&mut self) -> Result<()> {
        let h = self.step_size();
        self.calculate_state_k1();
        self.calculate_state_k2(h)?;
        self.calculate_state_k3(h)?;
        self.calculate_state_k4(h)?;
        self.rk45_state.calculate_add_terms();
        self.calculate_next_state(h)?;
        self.state = self.next.clone();

        self.points.push(Point {
            t: self.state.t,
            theta: self.state.theta,
            psip: self.state.psip,
            rho: self.state.rho,
            zeta: self.state.zeta,
        });

        Ok(())
    }

    fn step_size(&self) -> f64 {
        1e-1 // Hardcoded for now.
    }

    fn calculate_state_k1(&mut self) {
        self.rk45_state.k1[0] = self.state.theta_dot;
        self.rk45_state.k1[1] = self.state.psip_dot;
        self.rk45_state.k1[2] = self.state.rho_dot;
        self.rk45_state.k1[3] = self.state.zeta_dot;
    }

    fn calculate_state_k2(&mut self, h: f64) -> Result<()> {
        self.state2.t += h / 2.0;
        self.state2.theta += h / 2.0 * self.rk45_state.k1[0];
        self.state2.psip += h / 2.0 * self.rk45_state.k1[1];
        self.state2.rho += h / 2.0 * self.rk45_state.k1[2];
        self.state2.zeta += h / 2.0 * self.rk45_state.k1[3];
        self.state2
            .evaluate(self.qfactor, self.current, self.bfield)?;
        self.rk45_state.k2[0] = self.state2.theta_dot;
        self.rk45_state.k2[1] = self.state2.psip_dot;
        self.rk45_state.k2[2] = self.state2.rho_dot;
        self.rk45_state.k2[3] = self.state2.zeta_dot;

        Ok(())
    }

    fn calculate_state_k3(&mut self, h: f64) -> Result<()> {
        self.state3.t += h / 2.0;
        self.state3.theta += h / 2.0 * self.rk45_state.k2[0];
        self.state3.psip += h / 2.0 * self.rk45_state.k2[1];
        self.state3.rho += h / 2.0 * self.rk45_state.k2[2];
        self.state3.zeta += h / 2.0 * self.rk45_state.k2[3];
        self.state3
            .evaluate(self.qfactor, self.current, self.bfield)?;
        self.rk45_state.k3[0] = self.state3.theta_dot;
        self.rk45_state.k3[1] = self.state3.psip_dot;
        self.rk45_state.k3[2] = self.state3.rho_dot;
        self.rk45_state.k3[3] = self.state3.zeta_dot;
        Ok(())
    }

    fn calculate_state_k4(&mut self, h: f64) -> Result<()> {
        self.state4.t += h;
        self.state4.theta += h * self.rk45_state.k3[0];
        self.state4.psip += h * self.rk45_state.k3[1];
        self.state4.rho += h * self.rk45_state.k3[2];
        self.state4.zeta += h * self.rk45_state.k3[3];
        self.state4
            .evaluate(self.qfactor, self.current, self.bfield)?;
        self.rk45_state.k4[0] = self.state4.theta_dot;
        self.rk45_state.k4[1] = self.state4.psip_dot;
        self.rk45_state.k4[2] = self.state4.rho_dot;
        self.rk45_state.k4[3] = self.state4.zeta_dot;
        Ok(())
    }

    fn calculate_next_state(&mut self, h: f64) -> Result<()> {
        self.next.t += h;
        self.next.theta += h / 6.0 * self.rk45_state.k_add[0];
        self.next.psip += h / 6.0 * self.rk45_state.k_add[1];
        self.next.rho += h / 6.0 * self.rk45_state.k_add[2];
        self.next.zeta += h / 6.0 * self.rk45_state.k_add[3];
        self.next.evaluate(self.qfactor, self.current, self.bfield)
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use tokamak_equilibria::{Bfield, Current, Qfactor};

    use crate::*;

    #[test]
    fn test_system_init() {
        let path = PathBuf::from("./data.nc");
        let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
        let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
        let current = Current::from_dataset(&path, "akima").unwrap();
        // let initial = InitialConditions::new(0.1, 0.05, 0.01, 0.1);
        let initial = InitialConditions {
            t0: 0.0,
            theta0: 0.0,
            psip0: 0.01,
            rho0: 0.00,
            zeta0: 0.1,
            mu: 1e-6,
            pzeta: -0.01,
        };

        let mut system = System::init(&qfactor, &current, &bfield, initial).unwrap();

        for i in 0..4 {
            println!("step {i}, {}", system.state);
            match system.step() {
                Ok(_) => (),
                Err(err) => {
                    println!("{err:?}");
                    println!("Final State:");
                    println!("step {i}, {}", system.state);
                }
            }
        }
    }
}
