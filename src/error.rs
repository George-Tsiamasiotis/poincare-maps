#[derive(thiserror::Error, Debug)]
pub enum MapError {
    /// Error Initializing System.
    #[error("Error Initializing System")]
    SystemInitError,

    /// /// Equilibrium error.
    /// #[error("Equilibrium Error: {msg}")]
    /// EquilibriumError {
    ///     #[source]
    ///     source: tokamak_equilibria::EqError,
    ///     msg: Box<str>,
    /// },

    #[error("Equilibrium Error")]
    EquilibriumError(#[from] tokamak_equilibria::EqError), // TODO: fix msg

    /// Interpolaton Error.
    #[error("Interpolation Error: {msg}")]
    InterpolationError {
        #[source]
        source: rsl_interpolation::InterpolationError,
        msg: Box<str>,
    },
}
