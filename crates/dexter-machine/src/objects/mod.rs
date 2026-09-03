//! Machine representation objects.

pub(crate) mod getters;
pub(crate) mod nc_flux;

pub(crate) mod bfields;
pub(crate) mod currents;
pub(crate) mod flute_mode;
pub(crate) mod geometries;
pub(crate) mod nc_flute_mode;
pub(crate) mod perturbation;
pub(crate) mod qfactors;

/// Describes the type of machine.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MachineType {
    /// Machine reconstructed from numerical data.
    ///
    /// Evaluations are calculated by interpolation.
    Numerical,
    /// Machine described by analytical formulas.
    ///
    /// Evaluations are calculated by simply evaluating the formulas.
    Analytical,
}

/// Helper struct to define the Last Closed Flux Surface (LCFS) with respect to one of the
/// two fluxes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LastClosedFluxSurface {
    /// Last closed toroidal flux surface.
    Toroidal(f64),
    /// Last closed poloidal flux surface.
    Poloidal(f64),
}

impl LastClosedFluxSurface {
    /// Returns the LCFS value, regardless of its kind.
    #[must_use]
    pub fn value(&self) -> f64 {
        match *self {
            Self::Toroidal(value) | Self::Poloidal(value) => value,
        }
    }
}

/// Debug-asserts that all values of the slice are finite.
pub(crate) fn debug_assert_all_finite_values(values: &[f64]) {
    debug_assert!(
        values.iter().all(|element| element.is_finite()),
        "Non-finite values encountered."
    )
}

/// Debug-asserts that `r` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_r {
    ($r: expr) => {
        debug_assert!($r >= 0.0, "Encountered negative r")
    };
}

/// Debug-asserts that `psi` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_psi {
    ($psi: expr) => {
        debug_assert!($psi >= 0.0, "Encountered negative ψ")
    };
}

/// Debug-asserts that `psip` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_psip {
    ($psip: expr) => {
        debug_assert!($psip >= 0.0, "Encountered negative ψp")
    };
}

/// Debug-asserts that `theta` is 2π-modulo.
///
/// This is needed since the B array padding extends the interpolation domain.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_is_2pi_modulo {
    ($theta: expr) => {
        debug_assert!((0.0..=TAU).contains(&$theta), "Encountered non-2π-modulo θ")
    };
}

/// Debug-asserts that the wrapped expression `is_finite()` and returns it.
///
/// Use `std::debug_assert` and `std::dbg` syntax.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_is_finite {
    ($arg: expr) => {
        if std::cfg!(debug_assertions) {
            if $arg.is_finite() {
                $arg
            } else {
                panic!(
                    "[{}:{}:{}]: NaN value encountered",
                    std::file!(),
                    std::line!(),
                    std::column!(),
                )
            }
        } else {
            $arg
        }
    };
}
