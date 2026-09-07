//! Representation of the Trapped-Passing boundary curves on the `(E, Pζ, μ=const)` space.

use std::f64::consts::PI;

use dexter_machine::{FluxCoordinateState, Machine, MagneticFlux::*};
use ndarray::Array1;
use rsl_interpolation::{
    Accelerator, Accelerator2d, AkimaInterpolator, BuildInterpolator, Interpolation,
};

use crate::constants::TRAPPED_PASSING_BOUNDARY_DENSITY;

/// Representation of the Trapped-Passing boundary curves on the `(E, Pζ, μ=const)` space.
///
/// The Trapped-Passing boundary is defined by the curves:
///
/// 1. E = μB(ψp, 0)
/// 2. E = μB(ψp, π)
///
/// where `Pζ = -ψp`. Therefore, we can rewrite them in the `E-Pζ` plane.
///
/// 1. E(Pζ) = μB(-Pζ, 0)
/// 2. E(Pζ) = μB(-Pζ, π)
///
/// The peaks of the 2 wall parabolas and the axis parabola are always at the `Pζ/ψp_last = -1` and
/// `Pζ/ψp_last = 0` respectively. Therefore, the two curves lie in the `[-ψp_last, 0]` interval.
///
/// # Note
///
/// If the equilibrium has a `good` ψp coordinate, then it is used for all calculations, since
/// its faster. This also results in an equispaced [`TrappedPassingBoundary::pzeta_interval`].
///
/// If ψ is the only `good` coordinate, some extra conversions are necessary and the
/// `pzeta_interval` is no longer a linspace.
///
#[derive(Clone)]
pub(crate) struct TrappedPassingBoundary {
    /// The `Pζ = [-ψp_last, 0]` interval array.
    pub(crate) pzeta_interval: Box<[f64]>,
    /// The curve corresponding to the lower part of the boundary, defined by `θ=0` (eq. 1).
    pub(crate) lower: Box<[f64]>,
    /// The curve corresponding to the upper part of the boundary, defined by `θ=π` (eq. 2).
    pub(crate) upper: Box<[f64]>,
    /// The `Pζ` -> lower trapped-passing boundary interpolator.
    pub(crate) lower_interp: AkimaInterpolator,
    /// The `Pζ` -> upper trapped-passing boundary interpolator.
    pub(crate) upper_interp: AkimaInterpolator,
}

impl TrappedPassingBoundary {
    /// Calculates the two curves.
    #[must_use]
    pub(crate) fn new(machine: Machine, mu: f64) -> Self {
        // try `psip` first, since its faster.
        if machine.bfield().psip_state() == FluxCoordinateState::Good {
            Self::from_good_psip(machine, mu)
        } else if machine.bfield().psi_state() == FluxCoordinateState::Good {
            Self::from_good_psi(machine, mu)
        } else {
            unreachable!()
        }
    }

    /// Returns the [`TrappedPassingBoundary`]'s `Pζ = [-ψp_last, 0]` interval array.
    #[must_use]
    pub(crate) fn pzeta_interval(&self) -> Array1<f64> {
        Array1::from(self.pzeta_interval.clone())
    }

    /// Returns the [`TrappedPassingBoundary`]'s upper curve.
    #[must_use]
    pub(crate) fn upper(&self) -> Array1<f64> {
        Array1::from(self.upper.clone())
    }

    /// Returns the [`TrappedPassingBoundary`]'s lower curve.
    #[must_use]
    pub(crate) fn lower(&self) -> Array1<f64> {
        Array1::from(self.lower.clone())
    }
}

impl TrappedPassingBoundary {
    /// Returns `true` if an `(E, Pζ)` is inside the Trapped-Passing boundary.
    #[must_use]
    pub(crate) fn contains(&self, energy: f64, pzeta: f64, acc: &mut Accelerator) -> bool {
        !self.is_below(energy, pzeta, acc) && !self.is_above(energy, pzeta, acc)
    }

    /// Returns `true` if an `(E, Pζ)` is below the Trapped-Passing boundary's lower curve.
    #[must_use]
    pub(crate) fn is_below(&self, energy: f64, pzeta: f64, acc: &mut Accelerator) -> bool {
        let Some(lower_energy) = self
            .lower_interp
            .eval(&self.pzeta_interval, &self.lower, pzeta, acc)
            .ok()
        else {
            return false;
        };
        energy < lower_energy
    }

    /// Returns `true` if an `(E, Pζ)` is above the Trapped-Passing boundary's above curve.
    #[must_use]
    pub(crate) fn is_above(&self, energy: f64, pzeta: f64, acc: &mut Accelerator) -> bool {
        let Some(upper) = self
            .upper_interp
            .eval(&self.pzeta_interval, &self.upper, pzeta, acc)
            .ok()
        else {
            return false;
        };
        energy > upper
    }
}

impl TrappedPassingBoundary {
    /// Creates the boundary in a `Good ψp` equilibrium.
    ///
    /// This method cannot fail since `FluxState` is already checked and `ψp` is always in-bounds.
    fn from_good_psip(machine: Machine, mu: f64) -> Self {
        let psip_last = machine.qfactor().psip_last();
        let psip_interval =
            Array1::linspace(psip_last.value(), 0.0, TRAPPED_PASSING_BOUNDARY_DENSITY)
                .map(|value| Poloidal(*value));

        let acc = &mut Accelerator2d::new();

        let lower_array = mu
            * psip_interval.mapv(|psip| {
                machine
                    .bfield()
                    .eval_b(psip, 0.0, acc)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });
        let upper_array = mu
            * psip_interval.mapv(|psip| {
                machine
                    .bfield()
                    .eval_b(psip, PI, acc)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });

        let pzeta_interval = (psip_interval.map(|psip| -psip.value())).to_vec();
        let lower = lower_array.to_vec();
        let upper = upper_array.to_vec();
        let lower_interp = AkimaInterpolator::build(&pzeta_interval, &lower)
            .expect("sorted dataset and same shape by definition");
        let upper_interp = AkimaInterpolator::build(&pzeta_interval, &upper)
            .expect("sorted dataset and same shape by definition");

        Self {
            pzeta_interval: pzeta_interval.into_boxed_slice(),
            lower: lower.into_boxed_slice(),
            upper: upper.into_boxed_slice(),
            lower_interp,
            upper_interp,
        }
    }

    /// Creates the boundary in a `Good ψ` equilibrium.
    ///
    /// This method cannot fail since `FluxState` is already checked and `ψp/ψ` is always in-bounds.
    ///
    /// NOTE:
    /// We define the `psi_interval` first and the `psip_interval` second, since it is not
    /// guaranteed that `qfactor` defines `ψ(ψp)`. This results in a non-linspace `pzeta_interval`,
    /// but the values are correct.
    fn from_good_psi(machine: Machine, mu: f64) -> Self {
        let psi_last = machine.qfactor().psi_last();
        let psi_interval =
            Array1::linspace(psi_last.value(), 0.0, TRAPPED_PASSING_BOUNDARY_DENSITY)
                .map(|value| Toroidal(*value));

        let acc = &mut Accelerator2d::new();

        let lower_array = mu
            * psi_interval.mapv(|psi| {
                machine
                    .bfield()
                    .eval_b(psi, 0.0, acc)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });
        let upper_array = mu
            * psi_interval.mapv(|psi| {
                machine
                    .bfield()
                    .eval_b(psi, PI, acc)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });

        let psip_interval = psi_interval.mapv(|psi| {
            machine
                .qfactor()
                .eval_other(psi, acc.xacc())
                .expect("psi is always inbound and evaluation is defined")
        });

        let pzeta_interval = (psip_interval.map(|psip| -psip.value())).to_vec();
        let lower = lower_array.to_vec();
        let upper = upper_array.to_vec();
        let lower_interp = AkimaInterpolator::build(&pzeta_interval, &lower)
            .expect("sorted dataset and same shape by definition");
        let upper_interp = AkimaInterpolator::build(&pzeta_interval, &upper)
            .expect("sorted dataset and same shape by definition");

        Self {
            pzeta_interval: pzeta_interval.into_boxed_slice(),
            lower: lower.into_boxed_slice(),
            upper: upper.into_boxed_slice(),
            lower_interp,
            upper_interp,
        }
    }
}

impl std::fmt::Debug for TrappedPassingBoundary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TrappedPassingBoundary")
            .field("pzeta_interval", &self.pzeta_interval)
            .field("lower", &self.lower)
            .field("upper", &self.upper)
            .finish()
    }
}
