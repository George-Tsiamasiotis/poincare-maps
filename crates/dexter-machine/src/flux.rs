//! Defines the `MagneticFlux`.

use crate::EvalError;

/// Representation of the magnetic flux and its kind.
#[derive(Clone, Copy, PartialEq)]
pub enum MagneticFlux {
    /// Toroidal magnetic flux `ψ`.
    Toroidal(f64),
    /// Poloidal magnetic flux `ψp`.
    Poloidal(f64),
}

impl MagneticFlux {
    /// Returns the value of `self`, regardless the kind.
    pub fn value(&self) -> f64 {
        match *self {
            Self::Toroidal(value) | Self::Poloidal(value) => value,
        }
    }

    /// Returns the value of `self`, regardless the kind.
    pub fn value_mut(&mut self) -> &mut f64 {
        match *self {
            Self::Toroidal(ref mut value) => value,
            Self::Poloidal(ref mut value) => value,
        }
    }

    /// Returns the kind of `self` as a `Box<str>`.
    pub fn kind(&self) -> Box<str> {
        match self {
            Self::Toroidal(_) => "ψ".into(),
            Self::Poloidal(_) => "ψp".into(),
        }
    }

    /// Returns the flux's value if it is a [`MagneticFlux::Toroidal`] variant.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if `self` is a [`MagneticFlux::Poloidal`] variant.
    pub fn psi(&self) -> Result<f64, EvalError> {
        match self {
            Self::Toroidal(psi) => Ok(*psi),
            Self::Poloidal(_) => Err(EvalError::InvalidMagneticFlux),
        }
    }

    /// Returns the flux's value if it is a [`MagneticFlux::Poloidal`] variant.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if `self` is a [`MagneticFlux::Toroidal`] variant.
    pub fn psip(&self) -> Result<f64, EvalError> {
        match self {
            Self::Toroidal(_) => Err(EvalError::InvalidMagneticFlux),
            Self::Poloidal(psip) => Ok(*psip),
        }
    }
}

impl std::fmt::Debug for MagneticFlux {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::Toroidal(psi0) => write!(f, "ψ0: {psi0}"),
            Self::Poloidal(psip0) => write!(f, "ψp0: {psip0}"),
        }
    }
}

// ===============================================================================================

/// Compares discriminants and then forwards the comparison to the contained `f64`.
impl approx::AbsDiffEq for MagneticFlux {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        use std::mem::discriminant;
        if discriminant(self) != discriminant(other) {
            return false;
        }
        (self.value() - other.value()).abs() <= epsilon
    }
}

/// Compares discriminants and then forwards the comparison to the contained `f64`.
impl approx::RelativeEq for MagneticFlux {
    fn default_max_relative() -> Self::Epsilon {
        f64::EPSILON
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        use std::mem::discriminant;
        if discriminant(self) != discriminant(other) {
            return false;
        }
        self.value()
            .relative_eq(&other.value(), epsilon, max_relative)
    }
}

/// Compares discriminants and then forwards the comparison to the contained `f64`.
impl approx::UlpsEq for MagneticFlux {
    fn default_max_ulps() -> u32 {
        4
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        use std::mem::discriminant;
        if discriminant(self) != discriminant(other) {
            return false;
        }
        self.value().ulps_eq(&other.value(), epsilon, max_ulps)
    }
}

#[cfg(test)]
mod test {

    use crate::MagneticFlux::*;
    use approx::*;

    #[test]
    fn relative_eq() {
        assert_relative_eq!(Toroidal(0.01), Toroidal(0.01));
        assert_relative_eq!(Poloidal(0.01), Poloidal(0.01));

        assert_relative_ne!(Toroidal(0.01), Toroidal(0.02));
        assert_relative_ne!(Poloidal(0.01), Poloidal(0.02));
        assert_relative_ne!(Toroidal(0.01), Poloidal(0.01));
    }

    #[test]
    fn abs_diff_eq() {
        assert_abs_diff_eq!(Toroidal(0.01), Toroidal(0.01));
        assert_abs_diff_eq!(Poloidal(0.01), Poloidal(0.01));

        assert_abs_diff_ne!(Toroidal(0.01), Toroidal(0.02));
        assert_abs_diff_ne!(Poloidal(0.01), Poloidal(0.02));
        assert_abs_diff_ne!(Toroidal(0.01), Poloidal(0.01));
    }

    #[test]
    fn ulps_eq() {
        assert_ulps_eq!(Toroidal(0.01), Toroidal(0.01));
        assert_ulps_eq!(Poloidal(0.01), Poloidal(0.01));

        assert_ulps_ne!(Toroidal(0.01), Toroidal(0.02));
        assert_ulps_ne!(Poloidal(0.01), Poloidal(0.02));
        assert_ulps_ne!(Toroidal(0.01), Poloidal(0.01));
    }
}
