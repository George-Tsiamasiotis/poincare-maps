//! Representation of the COM space `(E, Pζ, μ=const)`.

#![expect(clippy::min_ident_chars, reason = "parabola a, b, c coefficients")]

use dexter_equilibrium::Equilibrium;
use ndarray::Array1;
use parabola::Parabola;
use rsl_interpolation::{Accelerator, Cache};

use crate::{COMError, COMs, TrappedPassingBoundary};

/// Representation of the COM space `(E, Pζ, μ=const)`.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct EnergyPzetaPlane {
    /// The magnetic axis (MA) parabola.
    axis_parabola: Parabola,
    /// The left wall (LW) parabola.
    left_wall_parabola: Parabola,
    /// The right wall (RW) parabola.
    right_wall_parabola: Parabola,
    /// The trapped-passing boundary curves.
    tp_boundary: TrappedPassingBoundary,
    /// The constant magnetic moment `μ`.
    mu: f64,
}

impl EnergyPzetaPlane {
    /// Creates a new `EnergyPzetaPlane` from a set of [`COMs`], in a given equilibrium.
    pub(crate) fn from_coms(equilibrium: &Equilibrium, coms: &COMs) -> Result<Self, COMError> {
        let Some(mu) = coms.mu else {
            return Err(COMError::UndefinedMu);
        };

        Ok(Self::from_mu(equilibrium, mu))
    }

    /// Creates a new `EnergyPzetaPlane` from a set magnetic moment `μ=const` value.
    pub(crate) fn from_mu(equilibrium: &Equilibrium, mu: f64) -> Self {
        Self {
            axis_parabola: Self::build_magnetic_axis_parabola(equilibrium, mu),
            left_wall_parabola: Self::build_left_wall_parabola(equilibrium, mu),
            right_wall_parabola: Self::build_right_wall_parabola(equilibrium, mu),
            tp_boundary: TrappedPassingBoundary::new(equilibrium, mu),
            mu,
        }
    }

    /// Constructs the Magnetic Axis (MA) parabola:
    ///
    /// E(Pζ) = Pζ²B²/2g² + μB.
    ///
    /// where all values are evaluated at `(ψ/ψp, θ) = (0, 0)`.
    #[must_use]
    fn build_magnetic_axis_parabola(equilibrium: &Equilibrium, mu: f64) -> Parabola {
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator::new();
        let cache = &mut Cache::<f64>::new();

        // Use `unwrap_or_else` for lazy evaluation.
        let gaxis = equilibrium.current.g_of_psi(0.0, acc1).unwrap_or_else(|_| {
            equilibrium
                .current
                .g_of_psip(0.0, acc1)
                .expect("At least one of the evaluations will always succeed")
        });
        let baxis = equilibrium
            .bfield // This might be redundant
            .b_of_psi(0.0, 0.0, acc1, acc2, cache)
            .unwrap_or_else(|_| {
                equilibrium
                    .bfield
                    .b_of_psip(0.0, 0.0, acc1, acc2, cache)
                    .expect("At least one of the evaluations will always succeed")
            });

        Parabola {
            a: (baxis / gaxis).powi(2) / 2.0,
            b: 0.0,
            c: mu * baxis,
        }
    }

    /// Constructs the Left Wall (LW) parabola:
    ///
    /// E(Pζ) = (Pζ + ψp)²B²/2g² + μB.
    ///
    /// where all values are evaluated at `(ψ/ψp, θ) = (ψlast/ψplast, π)`.
    #[must_use]
    fn build_left_wall_parabola(equilibrium: &Equilibrium, mu: f64) -> Parabola {
        use std::f64::consts::PI;
        let psi_last = equilibrium.qfactor.psi_last();
        let psip_last = equilibrium.qfactor.psip_last();
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator::new();
        let cache = &mut Cache::<f64>::new();

        // Use `unwrap_or_else` for lazy evaluation.
        let glast = equilibrium
            .current
            .g_of_psi(psi_last, acc1)
            .unwrap_or_else(|_| {
                equilibrium
                    .current
                    .g_of_psip(psip_last, acc1)
                    .expect("At least one of the evaluations will always succeed")
            });
        let blast = equilibrium
            .bfield
            .b_of_psi(psi_last, PI, acc1, acc2, cache)
            .unwrap_or_else(|_| {
                equilibrium
                    .bfield
                    .b_of_psip(psip_last, PI, acc1, acc2, cache)
                    .expect("At least one of the evaluations will always succeed")
            });

        let a = (blast / glast).powi(2) / 2.0;
        let h = psip_last;
        let k = mu * blast;
        Parabola::from_square(a, h, k)
    }

    /// Constructs the Right Wall (RW) parabola:
    ///
    /// E(Pζ) = (Pζ + ψp)²B²/2g² + μB.
    ///
    /// where all values are evaluated at `(ψ/ψp, θ) = (ψlast/ψplast, 0)`.
    #[must_use]
    fn build_right_wall_parabola(equilibrium: &Equilibrium, mu: f64) -> Parabola {
        let psi_last = equilibrium.qfactor.psi_last();
        let psip_last = equilibrium.qfactor.psip_last();
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator::new();
        let cache = &mut Cache::<f64>::new();

        // Use `unwrap_or_else` for lazy evaluation.
        let glast = equilibrium
            .current
            .g_of_psi(psi_last, acc1)
            .unwrap_or_else(|_| {
                equilibrium
                    .current
                    .g_of_psip(psip_last, acc1)
                    .expect("At least one of the evaluations will always succeed")
            });
        let blast = equilibrium
            .bfield
            .b_of_psi(psi_last, 0.0, acc1, acc2, cache)
            .unwrap_or_else(|_| {
                equilibrium
                    .bfield
                    .b_of_psip(psip_last, 0.0, acc1, acc2, cache)
                    .expect("At least one of the evaluations will always succeed")
            });

        let a = (blast / glast).powi(2) / 2.0;
        let h = psip_last;
        let k = mu * blast;
        Parabola::from_square(a, h, k)
    }
}

impl EnergyPzetaPlane {
    /// Returns a reference to the magnetic axis (MA) [`Parabola`].
    #[must_use]
    pub fn axis_parabola(&self) -> &Parabola {
        &self.axis_parabola
    }

    /// Returns a reference to the left wall (LW) [`Parabola`].
    #[must_use]
    pub fn left_wall_parabola(&self) -> &Parabola {
        &self.left_wall_parabola
    }

    /// Returns a reference to the right wall (RW) [`Parabola`].
    #[must_use]
    pub fn right_wall_parabola(&self) -> &Parabola {
        &self.right_wall_parabola
    }

    /// Returns a reference to the [`TrappedPassingBoundary`].
    #[must_use]
    pub fn tp_boundary(&self) -> &TrappedPassingBoundary {
        &self.tp_boundary
    }

    /// Returns the constant magnetic moment `μ`.
    #[must_use]
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Returns the [`TrappedPassingBoundary`]'s `Pζ = [-ψp_last, 0]` interval array.
    #[must_use]
    pub fn tp_pzeta_interval(&self) -> Array1<f64> {
        self.tp_boundary.pzeta_interval()
    }

    /// Returns the [`TrappedPassingBoundary`]'s upper curve.
    #[must_use]
    pub fn tp_upper(&self) -> Array1<f64> {
        self.tp_boundary.upper()
    }

    /// Returns the [`TrappedPassingBoundary`]'s lower curve.
    #[must_use]
    pub fn tp_lower(&self) -> Array1<f64> {
        self.tp_boundary.lower()
    }
}
