//! Implementation of the guiding center Hamiltonian and canonical equations.

#![allow(clippy::missing_docs_in_private_items, reason = "self-explanatory")]

use std::f64::consts::TAU;

use dexter_equilibrium::{Equilibrium, EvalError};

use crate::particle::IntegrationCaches;
use crate::{FluxCoordinate, InitialConditions, InitialFlux, SimulationError};

/// State of the Guiding Center at each step.
///
/// Stores all the intermediate values needed for the calculation of the final time derivatives.
///
/// Corresponds to a single specific point in configuration space, e.g. all values are calculated
/// at the same (t, ψ, θ, ζ, ρ, μ), or (t, ψp, θ, ζ, ρ, μ) point, depending on the value of the
/// `coordinate` field.
#[derive(Debug, Clone)]
#[expect(clippy::min_ident_chars, reason = "symbols in the Hamiltonian")]
pub(crate) struct GCState {
    pub(crate) coordinate: FluxCoordinate,
    pub(crate) t: f64,

    // Dynamic variables
    pub(crate) psi: f64,
    pub(crate) psip: f64,
    pub(crate) theta: f64,
    pub(crate) zeta: f64,
    pub(crate) rho: f64,
    pub(crate) mu: f64,

    // Final time derivatives
    pub(crate) psi_dot: f64,
    pub(crate) psip_dot: f64,
    pub(crate) theta_dot: f64,
    pub(crate) zeta_dot: f64,
    pub(crate) rho_dot: f64,
    pub(crate) mu_dot: f64,

    // Modulos of angles, to avoid recalculation
    mod_theta: f64,
    mod_zeta: f64,

    pub(crate) ptheta: f64,
    pub(crate) pzeta: f64,
    pub(crate) energy: f64,

    // Equilibrium quantities
    b: f64,
    q: f64,
    g: f64,
    i: f64,
    p: f64,
    db_dflux: f64,
    db_dtheta: f64,
    db_dzeta: f64,
    dg_dflux: f64,
    di_dflux: f64,
    dp_dflux: f64,
    dp_dtheta: f64,
    dp_dzeta: f64,
    dp_dt: f64,

    delta: f64,
    ddelta_dtheta: f64,
    ddelta_dzeta: f64,

    /// Constants D, K, C, F, as they appear in the equations of motion.
    dterm: f64,
    kterm: f64,
    cterm: f64,
    fterm: f64,

    /// The intermediate value (μ + ρ^2Β).
    mu_par: f64,
    /// The intermediate value [dψ or dψp derivatives].
    dflux_brace: f64,
    /// The intermediate value [theta derivatives].
    dtheta_brace: f64,
    /// The intermediate value [zeta derivatives].
    dzeta_brace: f64,
    /// The intermediate value ρ*B^2/D.
    rho_bsquared_d: f64,
    /// The intermediate value g/D.
    g_over_d: f64,
    /// The intermediate value I/D.
    i_over_d: f64,
}

/// Creation and evaluation.
impl GCState {
    /// Creates a new `GCState` from a set of [`InitialConditions`] and evaluates it.
    pub(crate) fn new(
        initial: &InitialConditions,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<Self, SimulationError> {
        let (psi, psip): (f64, f64);
        let coordinate: FluxCoordinate;
        match initial.flux0 {
            InitialFlux::Toroidal(psi0) => {
                (psi, psip) = (psi0, f64::NAN);
                coordinate = FluxCoordinate::Toroidal;
            }
            InitialFlux::Poloidal(psip0) => {
                (psi, psip) = (f64::NAN, psip0);
                coordinate = FluxCoordinate::Poloidal;
            }
        }
        let Some(rho0) = initial.rho0 else {
            unreachable!("rho0 must be initialized at this point")
        };
        Self {
            t: initial.t0,
            psi,
            psip,
            rho: rho0,
            theta: initial.theta0,
            zeta: initial.zeta0,
            mu: initial.mu0,
            coordinate,
            ..Default::default()
        }
        .into_evaluated(equilibrium, caches)
    }

    /// Performs all evaluations and calculation of intermediate quantities and final time derivatives.
    pub(crate) fn evaluate(
        &mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<(), SimulationError> {
        // First do all the interpolations
        self.calculate_modulos();
        self.calculate_other_flux(equilibrium, caches)?;
        self.calculate_qfactor_quantities(equilibrium, caches)?;
        self.calculate_current_quantities(equilibrium, caches)?;
        self.calculate_bfield_quantities(equilibrium, caches)?;
        self.calculate_perturbation_quantities(equilibrium, caches)?;

        // Multiply with `q` where needed, depending on the `FluxCoordinate`
        self.adjust_for_flux();

        // Then the intermediate quantities that only depend on the interpolated quantities
        self.calculate_canonical_momenta();
        self.calculate_delta_terms();
        self.calculate_capitals();
        self.calculate_mu_par();
        self.calculate_braces();
        self.calculate_extras();
        self.calculate_energy();

        // And finally the derivatives
        self.calculate_flux_dots();
        self.calculate_theta_dot();
        self.calculate_zeta_dot();
        self.calculate_rho_dot();
        self.calculate_mu_dot();
        Ok(())
    }

    /// Returns the state evaluated, consuming self.
    pub(crate) fn into_evaluated(
        mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<Self, SimulationError> {
        self.evaluate(equilibrium, caches)?;
        Ok(self)
    }
}

/// Field calculations.
impl GCState {
    fn calculate_modulos(&mut self) {
        self.mod_theta = self.theta.rem_euclid(TAU);
        self.mod_zeta = self.zeta.rem_euclid(TAU);
    }

    /// Calculates the non-coordinate flux, if it is defined.
    fn calculate_other_flux(
        &mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<(), SimulationError> {
        let qfactor = equilibrium.qfactor.as_ref();
        if self.coordinate == FluxCoordinate::Toroidal {
            self.psip = match qfactor.psip_of_psi(self.psi, &mut caches.psi_acc) {
                Ok(psip) => psip,
                Err(EvalError::UndefinedEvaluation(..)) => f64::NAN,
                Err(err) => return Err(err.into()),
            }
        } else {
            self.psi = match qfactor.psi_of_psip(self.psip, &mut caches.psip_acc) {
                Ok(psi) => psi,
                Err(EvalError::UndefinedEvaluation(..)) => f64::NAN,
                Err(err) => return Err(err.into()),
            }
        }
        Ok(())
    }

    fn calculate_qfactor_quantities(
        &mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<(), SimulationError> {
        let qfactor = equilibrium.qfactor.as_ref();
        if self.coordinate == FluxCoordinate::Toroidal {
            self.q = qfactor.q_of_psi(self.psi, &mut caches.psi_acc)?;
        } else {
            self.q = qfactor.q_of_psip(self.psip, &mut caches.psip_acc)?;
        };
        Ok(())
    }

    fn calculate_current_quantities(
        &mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<(), SimulationError> {
        let current = equilibrium.current.as_ref();
        if self.coordinate == FluxCoordinate::Toroidal {
            self.g = current.g_of_psi(self.psi, &mut caches.psi_acc)?;
            self.i = current.i_of_psi(self.psi, &mut caches.psi_acc)?;
            self.dg_dflux = current.dg_dpsi(self.psi, &mut caches.psi_acc)?;
            self.di_dflux = current.di_dpsi(self.psi, &mut caches.psi_acc)?;
        } else {
            self.g = current.g_of_psip(self.psip, &mut caches.psip_acc)?;
            self.i = current.i_of_psip(self.psip, &mut caches.psip_acc)?;
            self.dg_dflux = current.dg_dpsip(self.psip, &mut caches.psip_acc)?;
            self.di_dflux = current.di_dpsip(self.psip, &mut caches.psip_acc)?;
        }
        Ok(())
    }

    #[rustfmt::skip]
    fn calculate_bfield_quantities(
        &mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<(), SimulationError>
    {
        let bfield = equilibrium.bfield.as_ref();
        if self.coordinate == FluxCoordinate::Toroidal {
            self.b          = bfield.b_of_psi           (self.psi, self.mod_theta, &mut caches.psi_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dflux   = bfield.db_dpsi            (self.psi, self.mod_theta, &mut caches.psi_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dtheta  = bfield.db_of_psi_dtheta   (self.psi, self.mod_theta, &mut caches.psi_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
        } else {
            self.b          = bfield.b_of_psip          (self.psip, self.mod_theta, &mut caches.psip_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dflux   = bfield.db_dpsip           (self.psip, self.mod_theta, &mut caches.psip_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dtheta  = bfield.db_of_psip_dtheta  (self.psip, self.mod_theta, &mut caches.psip_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
        }
        self.db_dzeta = 0.0; // Axisymmetric configuration for now
        Ok(())
    }

    #[rustfmt::skip]
    fn calculate_perturbation_quantities(
        &mut self,
        equilibrium: &Equilibrium,
        caches: &mut IntegrationCaches,
    ) -> Result<(), SimulationError>
    {
        let perturbation = &equilibrium.perturbation;
        if self.coordinate == FluxCoordinate::Toroidal {
            self.p          = perturbation.p_of_psi         (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dflux   = perturbation.dp_dpsi          (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dtheta  = perturbation.dp_of_psi_dtheta (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dzeta   = perturbation.dp_of_psi_dzeta  (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dt      = perturbation.dp_of_psi_dt     (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
        } else {
            self.p          = perturbation.p_of_psip        (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dflux   = perturbation.dp_dpsip         (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dtheta  = perturbation.dp_of_psip_dtheta(self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dzeta   = perturbation.dp_of_psip_dzeta (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
            self.dp_dt      = perturbation.dp_of_psip_dt    (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.mode_caches)?;
        }
        Ok(())
    }

    /// Since we are solving White's equations of motion, which are expressed through (ψp, θ, ζ, ρ),
    /// to integrate with respect to `ψ` we must apply the chain rule wherever derivatives with
    /// respect to `ψp` appear, in order to obtain the derivatives with respect to `ψ`. This results
    /// to a multiplication with `q`, for which we adjust here.
    ///
    /// Chain rule:
    ///     dF/dψp = dF/dψ * dψ/dψp = q * dF/dψ.
    ///
    /// Still, the final expression of `flux_dot` is actually `psip_dot`, since this adjustment
    /// simply calculates the correct derivatives at the current point in the configuration space.
    /// Therefore, using the chain rule, we must also multiply the final `flux_dot` with `q` to
    /// obtain `psi_dot`.
    fn adjust_for_flux(&mut self) {
        if self.coordinate == FluxCoordinate::Toroidal {
            self.dg_dflux *= self.q;
            self.di_dflux *= self.q;
            self.dp_dflux *= self.q;
            self.db_dflux *= self.q;
        };
    }

    fn calculate_canonical_momenta(&mut self) {
        self.ptheta = self.psi + self.rho * self.i;
        self.pzeta = self.rho * self.g - self.psip;
    }

    fn calculate_delta_terms(&mut self) {
        self.delta = 0.0;
        self.ddelta_dtheta = 0.0;
        self.ddelta_dzeta = 0.0;
    }

    /// Calculates the matrix coefficients denoted with capital letters that appear in the
    /// perturbed equations of motion.
    ///
    /// If `ψ` is the coordinate, we must apply the chainrule in all g, I, p derivatives (with
    /// respect to `ψp`), which results in multiplying them with q.
    fn calculate_capitals(&mut self) {
        self.fterm = self.q + (self.rho + self.p) * self.di_dflux + self.i * self.dp_dflux;
        self.fterm -= self.q * self.rho * self.ddelta_dtheta; // δ term

        self.cterm = -1.0 + (self.rho + self.p) * self.dg_dflux + self.g * self.dp_dflux;
        self.cterm -= self.q * self.rho * self.ddelta_dzeta; // δ term

        self.kterm = self.g * self.dp_dtheta - self.i * self.dp_dzeta; // same without δ

        self.dterm = self.g * self.fterm - self.i * self.cterm;
        self.dterm += self.delta * self.q * self.kterm; // δ term
    }

    /// Calculates (μ + ρ^2B).
    fn calculate_mu_par(&mut self) {
        self.mu_par = self.mu + self.rho.powi(2) * self.b;
    }

    /// Calculates the brackets:
    ///     - `[mu_par*dB_d<flux> + dΦ_d<flux>]`
    ///     - `[mu_par*dB_dtheta + dΦ_dtheta]`
    /// where Φ = 0.
    ///
    /// Depending on the flux coordinate, it multiplies by q-factor where needed.
    fn calculate_braces(&mut self) {
        self.dflux_brace = self.mu_par * self.db_dflux;
        self.dtheta_brace = self.mu_par * self.db_dtheta;
        self.dzeta_brace = self.mu_par * self.db_dzeta;
    }

    /// Calculates intermediate values that appear many times in the time derivatives.
    fn calculate_extras(&mut self) {
        self.rho_bsquared_d = self.rho * self.b.powi(2) / self.dterm;
        self.g_over_d = self.g / self.dterm;
        self.i_over_d = self.i / self.dterm;
    }

    fn calculate_energy(&mut self) {
        self.energy = self.energy()
    }

    /// See [`Self::adjust_for_flux`].
    fn calculate_flux_dots(&mut self) {
        let flux_dot = self.kterm * self.rho_bsquared_d - self.g_over_d * self.dtheta_brace
            + self.i_over_d * self.dzeta_brace;
        self.psi_dot = flux_dot * self.q;
        self.psip_dot = flux_dot;
    }

    fn calculate_theta_dot(&mut self) {
        self.theta_dot = -self.cterm * self.rho_bsquared_d + self.g_over_d * self.dflux_brace;
        self.theta_dot += -self.delta * self.q / self.dterm * self.dp_dt;
    }

    fn calculate_zeta_dot(&mut self) {
        self.zeta_dot = self.fterm * self.rho_bsquared_d - self.i_over_d * self.dflux_brace;
        self.zeta_dot += self.delta * self.q / self.dterm * (self.dtheta_brace + self.dp_dt);
    }

    fn calculate_rho_dot(&mut self) {
        self.rho_dot = self.cterm / self.dterm * self.dtheta_brace
            - self.kterm / self.dterm * self.dflux_brace
            - self.fterm / self.dterm * self.dzeta_brace
            - self.dp_dt;
    }

    fn calculate_mu_dot(&mut self) {
        self.mu_dot = 0.0;
    }

    /// Returns the Energy of the State.
    pub(crate) fn energy(&self) -> f64 {
        let parallel = self.parallel_energy();
        let perpendicular = self.perpendicular_energy();
        parallel + perpendicular
    }

    /// Returns the parallel energy of the State.
    pub(crate) fn parallel_energy(&self) -> f64 {
        // Use the ρ expression here, since the g^2 in the denominator causes numerical instability
        // on configurations with g=0.
        (self.rho * self.b).powi(2) / 2.0
    }

    /// Returns the perpendicular energy of the State.
    pub(crate) fn perpendicular_energy(&self) -> f64 {
        self.mu * self.b
    }
}

// ===============================================================================================

impl GCState {
    /// Returns a mutable reference to the flux coordinate.
    pub(crate) fn flux(&mut self) -> &mut f64 {
        match self.coordinate {
            FluxCoordinate::Toroidal => &mut self.psi,
            FluxCoordinate::Poloidal => &mut self.psip,
        }
    }

    /// Returns the array with the final time derivatives.
    pub(crate) fn dots(&self) -> [f64; 5] {
        [
            match self.coordinate {
                FluxCoordinate::Toroidal => self.psi_dot,
                FluxCoordinate::Poloidal => self.psip_dot,
            },
            self.theta_dot,
            self.zeta_dot,
            self.rho_dot,
            self.mu_dot,
        ]
    }
}

impl Default for GCState {
    /// Set all derived quantities to NaN and use the corresponding methods to set them up.
    ///
    /// If NaNs manage to propagate, it means there's something wrong.
    fn default() -> Self {
        Self {
            coordinate: Default::default(),
            t: f64::NAN,
            psi: f64::NAN,
            psip: f64::NAN,
            theta: f64::NAN,
            zeta: f64::NAN,
            rho: f64::NAN,
            mu: f64::NAN,
            psi_dot: f64::NAN,
            psip_dot: f64::NAN,
            theta_dot: f64::NAN,
            zeta_dot: f64::NAN,
            rho_dot: f64::NAN,
            mu_dot: f64::NAN,
            mod_theta: f64::NAN,
            mod_zeta: f64::NAN,
            ptheta: f64::NAN,
            pzeta: f64::NAN,
            energy: f64::NAN,
            b: f64::NAN,
            q: f64::NAN,
            g: f64::NAN,
            i: f64::NAN,
            p: f64::NAN,
            dg_dflux: f64::NAN,
            di_dflux: f64::NAN,
            db_dflux: f64::NAN,
            db_dtheta: f64::NAN,
            db_dzeta: f64::NAN,
            dp_dflux: f64::NAN,
            dp_dtheta: f64::NAN,
            dp_dzeta: f64::NAN,
            dp_dt: f64::NAN,
            delta: f64::NAN,
            ddelta_dtheta: f64::NAN,
            ddelta_dzeta: f64::NAN,
            dterm: f64::NAN,
            kterm: f64::NAN,
            cterm: f64::NAN,
            fterm: f64::NAN,
            mu_par: f64::NAN,
            dflux_brace: f64::NAN,
            dtheta_brace: f64::NAN,
            dzeta_brace: f64::NAN,
            rho_bsquared_d: f64::NAN,
            g_over_d: f64::NAN,
            i_over_d: f64::NAN,
        }
    }
}
