use core::f64;
use std::f64::consts::TAU;

use equilibrium::HarmonicCache;
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use rsl_interpolation::{Accelerator, Cache};

use crate::Result;
use crate::{Distance, InitialConditions, MagneticMoment};
use crate::{Flux, Radians, Time};

/// State of the System at each step.
///
/// Corresponds to a single specific point in configuration space, e.g. all values are calculated
/// at the same `θ`, `ψp`, `ρ`, `ζ` point.
#[derive(Clone)]
pub struct State {
    /// The `ψp` coordinate [`Accelerator`].
    pub xacc: Accelerator,
    /// The `θ` coordinate [`Accelerator`].
    pub yacc: Accelerator,
    /// The 2d Interpolation's [`Cache`].
    pub spline2d_cache: Cache<f64>,
    /// Holds a [`HarmonicCache`] for each harmonic of the [`Perturbation`], caching the expensive
    /// trigonometric and interpolation values.
    pub hcache: Vec<HarmonicCache>,

    /// The time of evaluation.
    pub time: Time,

    /// The `θ` angle.
    pub theta: Radians,
    /// The poloidal magnetic flux `ψp`.
    pub psip: Flux,
    /// The parallel gyro radius `ρ`.
    pub rho: Distance,
    /// The `ζ` angle.
    pub zeta: Radians,

    /// The magnetic moment.
    pub mu: MagneticMoment,
    /// The toroidal magnetic flux `ψ`.
    pub psi: Flux,
    /// The canonical momentum `Pθ`,
    pub ptheta: f64,
    /// The canonical momentum `Pζ`,
    pub pzeta: f64,

    /// The `θ` angle time derivative.
    pub theta_dot: f64,
    /// The magnetic flux `ψ` time derivative.
    pub psip_dot: f64,
    /// The parallel gyro radius `ρ` time derivative.
    pub rho_dot: f64,
    /// The `ζ` angle time derivative.
    pub zeta_dot: f64,

    /// The magnetic field strength.
    pub b: f64,
    /// The safety factor `q`.
    pub q: f64,
    /// The toroidal plasma current.
    pub g: f64,
    /// The poloidal plasma current.
    pub i: f64,
    /// The perturbation (sum of harmonics)
    pub p: f64,

    /// The magnetic field strength derivative with respect to `θ`.
    pub db_dtheta: f64,
    /// The magnetic field strength derivative with respect to `ζ`. Should always be 0, since we
    /// are dealing with axisymmetric equilibria.
    pub db_dzeta: f64,
    /// The magnetic field strength derivative with respect to `ψ`.
    pub db_dpsip: f64,
    /// The toroidal plasma current derivative with respect to `ψ`
    pub dg_dpsip: f64,
    /// The poloidal plasma current derivative with respect to `ψ`
    pub di_dpsip: f64,
    /// The perturbation's (sum of harmonics) derivative with respect to `ψp`.
    pub dp_dpsip: f64,
    /// The perturbation's (sum of harmonics) derivative with respect to `θ`.
    pub dp_dtheta: f64,
    /// The perturbation's (sum of harmonics) derivative with respect to `ζ`.
    pub dp_dzeta: f64,
    /// The perturbation's (sum of harmonics) derivative with respect to `t`.
    pub dp_dt: f64,

    /// The `D` coefficient.
    pub dterm: f64,
    /// The `K` coefficient.
    pub kterm: f64,
    /// The `C` coefficient.
    pub cterm: f64,
    /// The `F` coefficient.
    pub fterm: f64,

    /// The intermediate value .
    pub mu_par: f64,
    /// The intermediate value [psip derivatives].
    pub psip_brace: f64,
    /// The intermediate value [theta derivatives].
    pub theta_brace: f64,
    /// The intermediate value [zeta derivatives].
    pub zeta_brace: f64,
    /// The intermediate value ρ*B^2/D.
    pub rho_bsquared_d: f64,
    /// The intermediate value g/D.
    pub g_over_d: f64,
    /// The intermediate value i/D.
    pub i_over_d: f64,

    pub mod_theta: f64,
    pub mod_zeta: f64,
}

impl State {
    /// Creates a non-evaluated [`State`] from an initial conditions set.
    pub fn from_initial(initial: &InitialConditions) -> Self {
        Self {
            time: initial.time0,
            theta: initial.theta0,
            psip: initial.psip0,
            rho: initial.rho0,
            zeta: initial.zeta0,
            mu: initial.mu,
            ..Default::default()
        }
    }

    /// Returns the state evaluated, consuming self.
    pub fn into_evaluated(
        mut self,
        qfactor: &Qfactor,
        currents: &Currents,
        bfield: &Bfield,
        perturbation: &Perturbation,
    ) -> Result<Self> {
        self.evaluate(qfactor, currents, bfield, perturbation)?;
        Ok(self)
    }

    /// Evaluation all quantites derived by (θ, ψp, ρ, ζ, μ)
    pub fn evaluate(
        &mut self,
        qfactor: &Qfactor,
        currents: &Currents,
        bfield: &Bfield,
        perturbation: &Perturbation,
    ) -> Result<()> {
        // First do all the interpolations.
        self.calculate_modulos();
        self.calculate_qfactor_quantities(qfactor)?;
        self.calculate_currents_quantities(currents)?;
        self.calculate_bfield_quantities(bfield)?;
        self.calculate_perturbation(perturbation)?;

        // Then intermediate quantities that only depend on the already calculated interpolations.
        self.calculate_capitals();
        self.calculate_mu_par();
        self.calculate_braces();
        self.calculate_extras();
        self.calculate_canonical_momenta();

        // And finally the time derivatives.
        self.calculate_theta_dot();
        self.calculate_psip_dot();
        self.calculate_rho_dot();
        self.calculate_zeta_dot();
        Ok(())
    }

    fn calculate_modulos(&mut self) {
        self.mod_theta = self.theta % TAU;
        self.mod_zeta = self.zeta % TAU;
    }

    fn calculate_qfactor_quantities(&mut self, qfactor: &Qfactor) -> Result<()> {
        self.psi = qfactor.psi(self.psip, &mut self.xacc)?;
        self.q = qfactor.q(self.psip, &mut self.xacc)?;
        Ok(())
    }

    fn calculate_currents_quantities(&mut self, currents: &Currents) -> Result<()> {
        self.i = currents.i(self.psip, &mut self.xacc)?;
        self.g = currents.g(self.psip, &mut self.xacc)?;
        self.di_dpsip = currents.di_dpsip(self.psip, &mut self.xacc)?;
        self.dg_dpsip = currents.dg_dpsip(self.psip, &mut self.xacc)?;
        Ok(())
    }

    fn calculate_bfield_quantities(&mut self, bfield: &Bfield) -> Result<()> {
        self.b = bfield.b(
            self.psip,
            self.mod_theta,
            &mut self.xacc,
            &mut self.yacc,
            &mut self.spline2d_cache,
        )?;
        self.db_dtheta = bfield.db_dtheta(
            self.psip,
            self.mod_theta,
            &mut self.xacc,
            &mut self.yacc,
            &mut self.spline2d_cache,
        )?;
        self.db_dpsip = bfield.db_dpsip(
            self.psip,
            self.mod_theta,
            &mut self.xacc,
            &mut self.yacc,
            &mut self.spline2d_cache,
        )?;
        self.db_dzeta = 0.0; // Axisymmetric configuration
        Ok(())
    }

    fn calculate_perturbation(&mut self, perturbation: &Perturbation) -> Result<()> {
        // This is necessary, since we can't know the number of harmonics from the start. The
        // hcache vec should be cloned in each solver state.
        if self.hcache.is_empty() {
            self.hcache = vec![HarmonicCache::default(); perturbation.harmonics.len()];
        }

        self.p = perturbation.p(
            self.psip,
            self.mod_theta,
            self.mod_zeta,
            &mut self.hcache,
            &mut self.xacc,
        )?;
        self.dp_dpsip = perturbation.dp_dpsip(
            self.psip,
            self.mod_theta,
            self.mod_zeta,
            &mut self.hcache,
            &mut self.xacc,
        )?;
        self.dp_dtheta = perturbation.dp_dtheta(
            self.psip,
            self.mod_theta,
            self.mod_zeta,
            &mut self.hcache,
            &mut self.xacc,
        )?;
        self.dp_dzeta = perturbation.dp_dzeta(
            self.psip,
            self.mod_theta,
            self.mod_zeta,
            &mut self.hcache,
            &mut self.xacc,
        )?;
        self.dp_dt = perturbation.dp_dt(
            self.psip,
            self.mod_theta,
            self.mod_zeta,
            &mut self.hcache,
            &mut self.xacc,
        )?;
        Ok(())
    }

    fn calculate_canonical_momenta(&mut self) {
        self.ptheta = self.psi + self.rho * self.i;
        self.pzeta = self.rho * self.g - self.psip;
    }

    /// Calculates the matrix coefficients denoted with capital letters that appear in the
    /// perturbed equations of motion.
    fn calculate_capitals(&mut self) {
        self.cterm = -1.0 + (self.rho + self.p) * self.dg_dpsip + self.g * self.dp_dpsip;
        self.fterm = self.q + (self.rho + self.p) * self.di_dpsip + self.i * self.dp_dpsip;
        self.kterm = self.g * self.dp_dtheta - self.i * self.dp_dzeta;
        self.dterm = self.g * self.fterm - self.i * self.cterm;
    }

    fn calculate_mu_par(&mut self) {
        self.mu_par = self.mu + self.rho.powi(2) * self.b;
    }

    fn calculate_braces(&mut self) {
        self.psip_brace = self.mu_par * self.db_dpsip;
        self.theta_brace = self.mu_par * self.db_dtheta;
        self.zeta_brace = self.mu_par * self.db_dzeta;
    }

    fn calculate_extras(&mut self) {
        self.rho_bsquared_d = self.rho * self.b.powi(2) / self.dterm;
        self.g_over_d = self.g / self.dterm;
        self.i_over_d = self.i / self.dterm;
    }

    fn calculate_theta_dot(&mut self) {
        self.theta_dot = -self.cterm * self.rho_bsquared_d + self.g_over_d * self.psip_brace;
    }

    fn calculate_psip_dot(&mut self) {
        self.psip_dot = self.kterm * self.rho_bsquared_d - self.g_over_d * self.theta_brace
            + self.i_over_d * self.zeta_brace;
    }

    fn calculate_rho_dot(&mut self) {
        self.rho_dot = self.cterm / self.dterm * self.theta_brace
            - self.kterm / self.dterm * self.psip_brace
            - self.fterm / self.dterm * self.zeta_brace
            - self.dp_dt;
    }

    fn calculate_zeta_dot(&mut self) {
        self.zeta_dot = self.fterm * self.rho_bsquared_d - self.i_over_d * self.psip_brace
    }

    pub fn energy(&self) -> f64 {
        let parallel = self.parallel_energy();
        let perpendicular = self.perpendicular_energy();
        parallel + perpendicular
    }

    pub fn parallel_energy(&self) -> f64 {
        (self.pzeta + self.psip).powi(2) * self.b.powi(2) / (2.0 * self.g.powi(2))
    }

    pub fn perpendicular_energy(&self) -> f64 {
        self.mu * self.b
    }
}

/// Helper struct for printing [`State`]'s independent variables.
pub struct Display {
    pub time: f64,
    pub theta: f64,
    pub psip: f64,
    pub rho: f64,
    pub zeta: f64,
    pub mu: f64,
}

impl std::fmt::Debug for Display {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("")
            .field("time", &format!("{:.5}", self.time))
            .field("theta", &format!("{:.5}", self.theta))
            .field("psip", &format!("{:.5}", self.psip))
            .field("rho", &format!("{:.5}", self.rho))
            .field("zeta", &format!("{:.5}", self.zeta))
            .field("mu", &format!("{:.5}", self.mu))
            .finish()
    }
}

impl Display {
    pub fn from_state(state: &State) -> Self {
        Self {
            time: state.time,
            theta: state.theta,
            psip: state.psip,
            rho: state.rho,
            zeta: state.zeta,
            mu: state.zeta,
        }
    }
}

impl Default for State {
    /// Set all derived quantities to NaN and use the corresponding methods to set them up
    fn default() -> Self {
        Self {
            xacc: Accelerator::new(),
            yacc: Accelerator::new(),
            spline2d_cache: Cache::new(),
            hcache: Vec::new(),
            time: f64::NAN,
            theta: f64::NAN,
            psip: f64::NAN,
            rho: f64::NAN,
            zeta: f64::NAN,
            mu: f64::NAN,
            psi: f64::NAN,
            ptheta: f64::NAN,
            pzeta: f64::NAN,
            theta_dot: f64::NAN,
            psip_dot: f64::NAN,
            rho_dot: f64::NAN,
            zeta_dot: f64::NAN,
            b: f64::NAN,
            q: f64::NAN,
            g: f64::NAN,
            i: f64::NAN,
            db_dtheta: f64::NAN,
            db_dzeta: f64::NAN,
            db_dpsip: f64::NAN,
            dg_dpsip: f64::NAN,
            di_dpsip: f64::NAN,
            dterm: f64::NAN,
            kterm: f64::NAN,
            cterm: f64::NAN,
            fterm: f64::NAN,
            p: f64::NAN,
            dp_dpsip: f64::NAN,
            dp_dtheta: f64::NAN,
            dp_dzeta: f64::NAN,
            dp_dt: f64::NAN,
            mu_par: f64::NAN,
            psip_brace: f64::NAN,
            theta_brace: f64::NAN,
            zeta_brace: f64::NAN,
            rho_bsquared_d: f64::NAN,
            g_over_d: f64::NAN,
            i_over_d: f64::NAN,
            mod_theta: f64::NAN,
            mod_zeta: f64::NAN,
        }
    }
}

/// Just to remove the Caches and Accelerators.
impl std::fmt::Debug for State {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("State")
            .field("time", &self.time)
            .field("theta", &self.theta)
            .field("psip", &self.psip)
            .field("rho", &self.rho)
            .field("zeta", &self.zeta)
            .field("mu", &self.mu)
            .field("theta_dot", &self.theta_dot)
            .field("psip_dot", &self.psip_dot)
            .field("rho_dot", &self.rho_dot)
            .field("zeta_dot", &self.zeta_dot)
            .field("b", &self.b)
            .field("q", &self.q)
            .field("g", &self.g)
            .field("i", &self.i)
            .field("p", &self.p)
            .field("db_dtheta", &self.db_dtheta)
            .field("db_dzeta", &self.db_dzeta)
            .field("db_dpsip", &self.db_dpsip)
            .field("dg_dpsip", &self.dg_dpsip)
            .field("di_dpsip", &self.di_dpsip)
            .field("dp_dpsip", &self.dp_dpsip)
            .field("dp_dtheta", &self.dp_dtheta)
            .field("dp_dzeta", &self.dp_dzeta)
            .field("dp_dt", &self.dp_dt)
            .field("dterm", &self.dterm)
            .field("kterm", &self.kterm)
            .field("cterm", &self.cterm)
            .field("fterm", &self.fterm)
            .field("mu_par", &self.mu_par)
            .field("psip_brace", &self.psip_brace)
            .field("theta_brace", &self.theta_brace)
            .field("zeta_brace", &self.zeta_brace)
            .field("rho_bsquared_d", &self.rho_bsquared_d)
            .field("g_over_d", &self.g_over_d)
            .field("i_over_d", &self.i_over_d)
            .finish()
    }
}
