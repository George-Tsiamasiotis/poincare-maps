pub(crate) struct Rk45State {
    pub(crate) k1: [f64; 4],
    pub(crate) k2: [f64; 4],
    pub(crate) k3: [f64; 4],
    pub(crate) k4: [f64; 4],
    pub(crate) k_add: [f64; 4],
}

impl Rk45State {
    pub(crate) fn calculate_add_terms(&mut self) {
        for i in 0..self.k_add.len() {
            self.k_add[i] = self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i];
        }
    }
}

impl Default for Rk45State {
    fn default() -> Self {
        Self {
            k1: [f64::NAN, f64::NAN, f64::NAN, f64::NAN],
            k2: [f64::NAN, f64::NAN, f64::NAN, f64::NAN],
            k3: [f64::NAN, f64::NAN, f64::NAN, f64::NAN],
            k4: [f64::NAN, f64::NAN, f64::NAN, f64::NAN],
            k_add: [f64::NAN, f64::NAN, f64::NAN, f64::NAN],
        }
    }
}
