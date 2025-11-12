//! The names each variable is expected to appear in the netCDF file.
//!
//! If the naming convention changes, this is the only file we must update.

// ================== Scalars ==================

/// Magnetic field strength on the axis `B0` **in \[T\]**.
pub const BAXIS: &str = "baxis";
/// The tokamak's major radius `R` **in \[m\]**.
pub const RAXIS: &str = "raxis";

// ================= Coordinates =================

/// The boozer toroidal angle `θ` **in \[rads\]**.
pub const THETA: &str = "theta";
/// The poloidal flux `ψp` **in Normalized Units**.
pub const PSIP: &str = "psip";
/// The toroidal flux `ψ` **in Normalized Units**.
pub const PSI: &str = "psi";
/// The poloidal mode numbers `m`.
pub const M: &str = "m";
/// The toroidal mode numbers `n`.
pub const N: &str = "n";

// ================ 1D Variables ================

/// q(ψp): The safety factor.
pub const Q: &str = "q";
/// g(ψp): The covariant toroidal plasma current **in \[Tm\]**.
pub const G: &str = "g";
/// I(ψp): The covariant poloidal plasma current **in \[Tm\]**.
pub const I: &str = "i";
/// g(ψp): The covariant toroidal plasma current **in Normalized Units**.
pub const G_NORM: &str = "g_norm";
/// I(ψp): The covariant poloidal plasma current **in Normalized Units**.
pub const I_NORM: &str = "i_norm";

// ================ 2D Variables ================

/// B(ψp, θ): The magnetic field strength in **in \[T\]**.
pub const B: &str = "b";
/// B(ψp, θ): The magnetic field strength in **in Normalized Units**.
pub const B_NORM: &str = "b_norm";
/// R(ψp, θ): The `R` coordinate with respect to boozer coordinates **in \[m\]**.
pub const R: &str = "R";
/// Z(ψp, θ): The `Z` coordinate with respect to boozer coordinates **in \[m\]**.
pub const Z: &str = "Z";

// ================ 3D Variables ================

/// The 3D array containing all the `α{m,n}(ψp)` 1D arrays **in Normalized Units**.
pub const ALPHAS_NORM: &str = "alphas_norm";
/// The 3D array containing all the `α{m,n}(ψp)` 1D arrays **i \[m\]**.
pub const ALPHAS: &str = "alphas";
/// The 3D array containing all the `φ{m,n}(ψp)` 1D arrays.
pub const PHASES: &str = "phases";
