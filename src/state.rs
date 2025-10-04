use core::f64;

use crate::{Bfield, Current, Qfactor};
use rsl_interpolation::Accelerator;

use crate::{InitialConditions, Result};

/// State of the System at each step.
///
/// Corresponds to a single specific point in configuration space, e.g. all values are calculated
/// at the same `θ`, `ψ_p`, `ρ`, `ζ` point.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct State {
    /// The `ψ_p` coordinate [`Accelerator`].
    xacc: Accelerator,
    /// The `θ_p` coordinate [`Accelerator`].
    yacc: Accelerator,

    /// The time of evaluation.
    pub(crate) t: f64,

    /// The `θ` angle.
    pub(crate) theta: f64,
    /// The poloidal magnetic flux `ψ_p`.
    pub(crate) psip: f64,
    /// The parallel gyro radius `ρ`.
    pub(crate) rho: f64,
    /// The `ζ` angle.
    pub(crate) zeta: f64,

    /// The magnetic moment.
    pub(crate) mu: f64,
    /// The toroidal magnetic flux `ψ`.
    pub(crate) psi: f64,
    /// The canonical momentum `Pθ`,
    pub(crate) ptheta: f64,
    /// The canonical momentum `Pζ`,
    pub(crate) pzeta: f64,
    /// Not sure yet, associated with the perturbations.
    pub(crate) a: f64,

    /// The `θ` angle time derivative.
    pub(crate) theta_dot: f64,
    /// The magnetic flux `ψ` time derivative.
    pub(crate) psip_dot: f64,
    /// The parallel gyro radius `ρ` time derivative.
    pub(crate) rho_dot: f64,
    /// The `ζ` angle time derivative.
    pub(crate) zeta_dot: f64,

    /// The magnetic field strength.
    pub(crate) b: f64,
    /// The safety factor `q`.
    pub(crate) q: f64,
    /// The toroidal plasma current.
    pub(crate) g: f64,
    /// The poloidal plasma current.
    pub(crate) i: f64,

    /// The magnetic field strength derivative with respect to `θ`.
    pub(crate) db_dtheta: f64,
    /// The magnetic field strength derivative with respect to `ζ`. Should always be 0, since we
    /// are dealing with axisymmetric equilibria.
    pub(crate) db_dzeta: f64,
    /// The magnetic field strength derivative with respect to `ψ`.
    pub(crate) db_dpsip: f64,
    /// The toroidal plasma current derivative with respect to `ψ`
    pub(crate) dg_dpsip: f64,
    /// The poloidal plasma current derivative with respect to `ψ`
    pub(crate) di_dpsip: f64,
    /// Not sure yet, associated with the perturbations.
    pub(crate) da_dpsip: f64,
    /// Not sure yet, associated with the perturbations.
    pub(crate) da_dtheta: f64,
    /// Not sure yet, associated with the perturbations.
    pub(crate) da_dzeta: f64,

    /// The `D` coefficient.
    pub(crate) dterm: f64,
    /// The `K` coefficient.
    pub(crate) kterm: f64,
    /// The `C` coefficient.
    pub(crate) cterm: f64,
    /// The `F` coefficient.
    pub(crate) fterm: f64,

    /// The intermediate value (μ+ρ^2Β).
    pub(crate) mu_par: f64,
    /// The intermediate value [<psip derivatives>].
    pub(crate) psip_brace: f64,
    /// The intermediate value [<theta derivatives>].
    pub(crate) theta_brace: f64,
    /// The intermediate value [<zeta derivatives>].
    pub(crate) zeta_brace: f64,
    /// The intermediate value ρ*B^2/D.
    pub(crate) rho_bsquared_d: f64,
    /// The intermediate value g/D.
    pub(crate) g_over_d: f64,
    /// The intermediate value i/D.
    pub(crate) i_over_d: f64,
}

impl State {
    pub(crate) fn init(initial: &InitialConditions) -> Result<Self> {
        // Set all derived quantities to NaN and use the corresponding methods to set them up
        Ok(Self {
            xacc: Accelerator::new(),
            yacc: Accelerator::new(),
            t: initial.t0,
            theta: initial.theta0,
            psip: initial.psip0,
            rho: initial.rho0,
            zeta: initial.zeta0,
            mu: initial.mu,
            psi: f64::NAN,
            ptheta: f64::NAN,
            pzeta: initial.pzeta,
            theta_dot: f64::NAN,
            psip_dot: f64::NAN,
            rho_dot: f64::NAN,
            zeta_dot: f64::NAN,
            b: f64::NAN,
            q: f64::NAN,
            g: f64::NAN,
            i: f64::NAN,
            db_dtheta: f64::NAN,
            db_dzeta: 0.0,
            db_dpsip: f64::NAN,
            dg_dpsip: f64::NAN,
            di_dpsip: f64::NAN,
            dterm: f64::NAN,
            kterm: f64::NAN,
            cterm: f64::NAN,
            fterm: f64::NAN,
            a: f64::NAN,
            da_dpsip: f64::NAN,
            da_dtheta: f64::NAN,
            da_dzeta: f64::NAN,
            mu_par: f64::NAN,
            psip_brace: f64::NAN,
            theta_brace: f64::NAN,
            zeta_brace: f64::NAN,
            rho_bsquared_d: f64::NAN,
            g_over_d: f64::NAN,
            i_over_d: f64::NAN,
        })
    }

    #[allow(dead_code)]
    pub(crate) fn evaluate(
        &mut self,
        qfactor: &Qfactor,
        current: &Current,
        bfield: &Bfield,
    ) -> Result<()> {
        // First do all the interpolations.
        self.calculate_qfactor_quantities(qfactor)?;
        self.calculate_current_quantities(current)?;
        self.calculate_bfield_quantities(bfield)?;

        self.calculate_perturbation()?; // TODO: when?

        // Then intermediate quantities that only depend on the already calculated interpolations.
        self.calculate_capitals();
        self.calculate_mu_par();
        self.calculate_braces();
        self.calculate_extras();

        self.calculate_theta_dot();
        self.calculate_psip_dot();
        self.calculate_rho_dot();
        self.calculate_zeta_dot();
        Ok(())
    }

    fn calculate_qfactor_quantities(&mut self, qfactor: &Qfactor) -> Result<()> {
        self.q = qfactor.q(self.psip, &mut self.xacc)?;
        Ok(())
    }

    fn calculate_current_quantities(&mut self, current: &Current) -> Result<()> {
        self.i = current.i(self.psip, &mut self.xacc)?;
        self.g = current.g(self.psip, &mut self.xacc)?;
        self.di_dpsip = current.di_dpsip(self.psip, &mut self.xacc)?;
        self.dg_dpsip = current.dg_dpsip(self.psip, &mut self.xacc)?;
        Ok(())
    }

    fn calculate_bfield_quantities(&mut self, bfield: &Bfield) -> Result<()> {
        self.b = bfield.b(self.psip, self.theta, &mut self.xacc, &mut self.yacc)?;
        self.db_dtheta = bfield.db_dtheta(self.psip, self.theta, &mut self.xacc, &mut self.yacc)?;
        self.db_dpsip = bfield.db_dpsip(self.psip, self.theta, &mut self.xacc, &mut self.yacc)?;
        Ok(())
    }

    fn calculate_perturbation(&mut self) -> Result<()> {
        self.a = 0.0;
        self.da_dpsip = 0.0;
        self.da_dtheta = 0.0;
        self.da_dzeta = 0.0;
        Ok(())
    }

    #[allow(dead_code)]
    fn calculate_canonical_momenta(&mut self) {
        self.ptheta = self.psi + self.rho * self.i; // TODO: find a way to calculate ψ
    }

    /// Calculates the matrix coefficients denoted with capital letters that appear in the
    /// perturbed equations of motion.
    fn calculate_capitals(&mut self) {
        self.cterm = -1.0 + (self.rho + self.a) * self.dg_dpsip + self.g * self.da_dpsip;
        self.fterm = self.q + (self.rho + self.a) * self.di_dpsip + self.i * self.da_dpsip;
        self.kterm = self.g * self.da_dtheta - self.i * self.da_dzeta;
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
        self.theta_dot = -self.cterm * self.rho_bsquared_d + self.g / self.dterm * self.psip_brace;
    }

    fn calculate_psip_dot(&mut self) {
        self.psip_dot = self.kterm * self.rho_bsquared_d - self.g_over_d * self.theta_brace
            + self.i_over_d * self.zeta_brace;
    }

    fn calculate_rho_dot(&mut self) {
        self.rho_dot = self.cterm / self.dterm * self.theta_brace
            - self.kterm / self.dterm * self.psip_brace
            - self.fterm / self.dterm * self.zeta_brace;
    }

    fn calculate_zeta_dot(&mut self) {
        self.zeta_dot = self.fterm * self.rho_bsquared_d - self.i_over_d * self.psip_brace
    }

    #[allow(dead_code)]
    pub(crate) fn energy(&self) -> f64 {
        let parallel = (self.pzeta + self.psip).powi(2) * self.b.powi(2) / (2.0 * self.g.powi(2));
        let perpandicular = self.mu * self.b;
        parallel + perpandicular
    }
}

impl std::fmt::Display for State {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("State Variables")
            .field("θ", &self.theta)
            .field("ψp", &self.psip)
            .field("ρ", &self.rho)
            .field("ζ", &self.zeta)
            .finish_non_exhaustive()
    }
}
