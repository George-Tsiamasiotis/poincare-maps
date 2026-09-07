//! Representation of the COM space `(E, Pζ, μ=const)`.

#![expect(unreachable_pub, reason = "api needs re-write")]
#![expect(clippy::min_ident_chars, reason = "parabola a, b, c coefficients")]

use dexter_machine::MagneticFlux::*;
use ndarray::Array1;
use parabola::Parabola;
use rsl_interpolation::Accelerator2d;

use dexter_machine::Machine;

use crate::COMError;
use crate::coms::{COMs, TrappedPassingBoundary};

/// Representation of the COM space `(E, Pζ, μ=const)`.
#[derive(Debug, Clone)]
pub(crate) struct EnergyPzetaPlane {
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
    pub(crate) fn from_coms(machine: Machine, coms: &COMs) -> Result<Self, COMError> {
        let Some(mu) = coms.mu else {
            return Err(COMError::UndefinedMu);
        };

        Ok(Self::from_mu(machine, mu))
    }

    /// Creates a new `EnergyPzetaPlane` from a set magnetic moment `μ=const` value.
    pub(crate) fn from_mu(objects: Machine, mu: f64) -> Self {
        Self {
            axis_parabola: Self::build_magnetic_axis_parabola(objects, mu),
            left_wall_parabola: Self::build_left_wall_parabola(objects, mu),
            right_wall_parabola: Self::build_right_wall_parabola(objects, mu),
            tp_boundary: TrappedPassingBoundary::new(objects, mu),
            mu,
        }
    }

    /// Constructs the Magnetic Axis (MA) parabola:
    ///
    /// E(Pζ) = Pζ²B²/2g² + μB.
    ///
    /// where all values are evaluated at `(ψ/ψp, θ) = (0, 0)`.
    #[must_use]
    fn build_magnetic_axis_parabola(objects: Machine, mu: f64) -> Parabola {
        let acc = &mut Accelerator2d::new();
        let psi_axis = Toroidal(0.0);
        let psip_axis = Poloidal(0.0);

        // Use `unwrap_or_else` for lazy evaluation.
        let gaxis = objects
            .current()
            .eval_g(psi_axis, acc.xacc())
            .unwrap_or_else(|_| {
                objects
                    .current()
                    .eval_g(psip_axis, acc.xacc())
                    .expect("At least one of the evaluations will always succeed")
            });
        let baxis = objects
            .bfield() // This might be redundant
            .eval_b(psi_axis, 0.0, acc)
            .unwrap_or_else(|_| {
                objects
                    .bfield()
                    .eval_b(psip_axis, 0.0, acc)
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
    fn build_left_wall_parabola(objects: Machine, mu: f64) -> Parabola {
        use std::f64::consts::PI;
        let psi_last = objects.qfactor().psi_last();
        let psip_last = objects.qfactor().psip_last();
        let acc = &mut Accelerator2d::new();

        // Use `unwrap_or_else` for lazy evaluation.
        let glast = objects
            .current()
            .eval_g(psi_last, acc.xacc())
            .unwrap_or_else(|_| {
                objects
                    .current()
                    .eval_g(psi_last, acc.xacc())
                    .expect("At least one of the evaluations will always succeed")
            });
        let blast = objects
            .bfield()
            .eval_b(psi_last, PI, acc)
            .unwrap_or_else(|_| {
                objects
                    .bfield()
                    .eval_b(psip_last, PI, acc)
                    .expect("At least one of the evaluations will always succeed")
            });

        let a = (blast / glast).powi(2) / 2.0;
        let h = psip_last.value();
        let k = mu * blast;
        Parabola::from_square(a, h, k)
    }

    /// Constructs the Right Wall (RW) parabola:
    ///
    /// E(Pζ) = (Pζ + ψp)²B²/2g² + μB.
    ///
    /// where all values are evaluated at `(ψ/ψp, θ) = (ψlast/ψplast, 0)`.
    #[must_use]
    fn build_right_wall_parabola(objects: Machine, mu: f64) -> Parabola {
        let psi_last = objects.qfactor().psi_last();
        let psip_last = objects.qfactor().psip_last();
        let acc = &mut Accelerator2d::new();

        // Use `unwrap_or_else` for lazy evaluation.
        let glast = objects
            .current()
            .eval_g(psi_last, acc.xacc())
            .unwrap_or_else(|_| {
                objects
                    .current()
                    .eval_g(psip_last, acc.xacc())
                    .expect("At least one of the evaluations will always succeed")
            });
        let blast = objects
            .bfield()
            .eval_b(psi_last, 0.0, acc)
            .unwrap_or_else(|_| {
                objects
                    .bfield()
                    .eval_b(psip_last, 0.0, acc)
                    .expect("At least one of the evaluations will always succeed")
            });

        let a = (blast / glast).powi(2) / 2.0;
        let h = psip_last.value();
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
