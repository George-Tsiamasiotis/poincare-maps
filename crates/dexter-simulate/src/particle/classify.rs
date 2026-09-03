//! Classification of a particle's orbit by projection on the `(E, Pζ, μ)` space.

use parabola::Point;
use rsl_interpolation::Accelerator;

use dexter_machine::Machine;

use crate::OrbitType;
use crate::coms::EnergyPzetaPlane;
use crate::particle::{IntegrationCaches, Particle};
use crate::state::GCState;

/// The position of an `(E-Pζ)` point on the `(E-Pζ)` plane, relative to the orbit
/// classification curves.
///
/// See the diagram for explanation.
#[non_exhaustive]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[expect(missing_docs, reason = "see diagram")]
pub enum EnergyPzetaPosition {
    /// No classification has been attempted.
    #[default]
    Undefined,
    Alpha,
    Beta,
    Gamma,
    Delta,
    Epsilon,
    Zeta,
    Eta,
    Theta,
    Iota,
    Kappa,
    Lambda,
    Mu,
    /// Not falling under any of the above categories.
    Unclassified,
}

impl EnergyPzetaPosition {
    /// Creates a new `EnergyPzetaPosition` from an `(E, Pζ)` plane and a point.
    #[expect(clippy::redundant_else, reason = "hopeless")]
    pub(crate) fn new(point: Point, plane: &EnergyPzetaPlane, theta0_dot: f64) -> Self {
        let psip_last = -plane
            .left_wall_parabola()
            .axis()
            .expect("parabola's 'a' checked");

        let pzeta = point.x;
        let energy = point.y;

        let mut acc = Accelerator::new();
        let is_in_axis = plane.axis_parabola().contains(point);
        let is_in_left_wall = plane.left_wall_parabola().contains(point);
        let is_in_right_wall = plane.right_wall_parabola().contains(point);
        let is_in_psip = (-psip_last..=0.0).contains(&pzeta);

        // Short-circuit condition to avoid the more expensive interpolations
        let might_be_trapped = !is_in_left_wall && is_in_psip;

        let tpb = plane.tp_boundary();
        let is_above_tp = might_be_trapped && tpb.is_above(energy, pzeta, &mut acc);
        let is_below_tp = might_be_trapped && tpb.is_below(energy, pzeta, &mut acc);
        let is_in_tpb = is_in_psip && !is_above_tp && !is_below_tp;

        // Inside left wall

        if is_in_left_wall {
            if is_in_axis {
                return Self::Unclassified; // No allowed orbits here.
            } else {
                return Self::Alpha;
            }
        }

        // Inside right wall, outside left wall. No need to check for left wall again

        if is_in_right_wall {
            if pzeta <= -psip_last {
                return Self::Beta;
            } else if is_in_tpb {
                return Self::Gamma; // Trapped-Lost
            } else if is_in_axis {
                return Self::Zeta;
            } else {
                if theta0_dot.is_sign_positive() {
                    return Self::Delta;
                } else {
                    return Self::Epsilon;
                }
            }
        }

        // Inside magnetic axis only. No need to check for walls again

        if is_in_axis {
            if is_in_tpb {
                return Self::Theta; // Potato
            } else {
                return Self::Eta;
            }
        }

        // Outside all parabolas. No need to check for any of them.

        if is_in_tpb {
            return Self::Iota; // Trapped-Confined
        }

        if is_above_tp {
            if theta0_dot.is_sign_positive() {
                return Self::Kappa;
            } else {
                return Self::Lambda;
            }
        } else {
            if pzeta > -psip_last {
                return Self::Mu; // Stagnated
            }
        }

        Self::Unclassified
    }

    /// Calculates the [`OrbitType`] from the resolved [`EnergyPzetaPosition`].
    #[expect(clippy::match_same_arms, reason = "clearer this way")]
    fn orbit_type(&self) -> OrbitType {
        match *self {
            Self::Undefined => OrbitType::Undefined,
            Self::Alpha => OrbitType::CuPassingConfined,
            Self::Beta => OrbitType::CuPassingLost,
            Self::Gamma => OrbitType::TrappedLost,
            Self::Delta => OrbitType::CoPassingLost,
            Self::Epsilon => OrbitType::CuPassingConfined,
            Self::Zeta => OrbitType::CoPassingLost,
            Self::Eta => OrbitType::CoPassingConfined,
            Self::Theta => OrbitType::Potato,
            Self::Iota => OrbitType::TrappedConfined,
            Self::Kappa => OrbitType::CoPassingConfined,
            Self::Lambda => OrbitType::CuPassingConfined,
            Self::Mu => OrbitType::Stagnated,
            Self::Unclassified => OrbitType::Unclassified,
        }
    }
}

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`OrbitType`] variant for each valid orbit calculation (including "erroneous" orbits),
/// and set [`OrbitType::Failed`] for hard errors.
///
/// If a `_plane` is passed, then that plane is used for the calculations, avoiding the
/// need for generating a new one. See [`Particle::_classify`].
pub(super) fn classify(
    particle: &mut Particle,
    machine: Machine,
    _plane: Option<&EnergyPzetaPlane>,
) {
    // =============== Particle Setup

    // Do not alter the `evolution` or `integration_status`
    let mut caches = IntegrationCaches {
        mode_caches: machine.perturbation().generate_caches(),
        ..Default::default()
    };

    // Return early if the initial flux happens to be exactly 0.0 or out of bounds.
    if particle.initial_conditions().flux0.value() == 0.0 {
        particle.orbit_type = OrbitType::Undefined;
        return;
    }
    if particle.initial_conditions.finalize(machine).is_err() {
        particle.orbit_type = OrbitType::Undefined;
        return;
    }
    let Ok(initial_state) = GCState::new(&particle.initial_conditions, machine, &mut caches) else {
        particle.orbit_type = OrbitType::Undefined;
        return;
    };

    particle.initial_energy = Some(initial_state.energy());
    particle.orbit_type = OrbitType::Unclassified; // Fallback

    // =============== Energy-Pζ plane Setup

    let mu = particle.initial_conditions.mu0;
    let plane = match _plane {
        Some(plane) => {
            #[expect(clippy::float_cmp, reason = "we need bit-to-bit equivalence")]
            if mu != plane.mu() {
                unreachable!("New EnergyPzetaPlanes must be generated");
            }
            plane
        }
        None => &EnergyPzetaPlane::from_mu(machine, mu),
    };

    check_parabola_alphas(plane);

    let pzeta = particle
        .initial_conditions
        .pzeta0
        .expect("Initial conditions have been finalized");
    let energy = particle
        .initial_energy
        .expect("initial energy has been calculated");
    let point = Point {
        x: pzeta,
        y: energy,
    };

    // =============== Routine

    particle.energy_pzeta_position =
        EnergyPzetaPosition::new(point, plane, initial_state.theta_dot);
    particle.orbit_type = particle.energy_pzeta_position.orbit_type();
}

// ===============================================================================================

/// Checks that all parabolas have nonzero `α` coefficients.
///
/// This is a corner case that should be fatal.
fn check_parabola_alphas(plane: &EnergyPzetaPlane) {
    assert!(
        plane.axis_parabola().a != 0.0,
        "Encountered zero 'a' coefficient in magnetic axis parabola"
    );
    assert!(
        plane.left_wall_parabola().a != 0.0,
        "Encountered zero 'a' coefficient in left wall parabola"
    );
    assert!(
        plane.right_wall_parabola().a != 0.0,
        "Encountered zero 'a' coefficient in right wall parabola"
    );
}
