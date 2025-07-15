//! Wrapper functions around [`rust-GSL`]'s wrapper functions. This crate
//! provides a cleaner interface, but also enables different splines to use the same [`Accelerator`].
//! This provides a significant performance boost when evalulating many splines with the same
//! datapoints at the same point, a case which comes up a lot in many physics calculations, such as
//! ODE problems.
//!
//! ## Example
//! ```
//! # use gsl_splines::{Accelerator, InterpolationType, Spline, SplineError};
//! # use ndarray::Array1;
//! #
//! # fn main() -> Result<(), SplineError>{
//! // Data creation
//! let xdata = Array1::linspace(0.0, 3.0, 100);
//! let ydata1 = xdata.sin();
//! let ydata2 = xdata.cos();
//!
//! // Interpolation type and Accelerator
//! let typ = InterpolationType::Cubic;
//! let acc = Accelerator::new();
//!
//! // Spline creation
//! let mut spline1 = Spline::build(typ, &xdata, &ydata1, acc.clone())?;
//! let mut spline2 = Spline::build(typ, &xdata, &ydata2, acc.clone())?;
//!
//! // 2 Different splines evaluating on a different point. The Accelerator's cache moves back and
//! // forth. In this case it is better to use 2 seperate Accelerators.
//! for x in Array1::linspace(1.0, 1.0001, 9) {
//!     spline1.eval(x)?;
//!     spline2.eval(x + 1.0)?;
//! }
//!
//! println!("{:?}", spline1.acc);
//! println!("{:?}", spline2.acc);
//!
//! spline1.reset_acc();
//! spline2.reset_acc();
//! println!();
//!
//! // 2 Different splines evaluating on the same point. Both splines can use the same Accelerator,
//! // so the index is searched only once.
//! for x in Array1::linspace(1.0, 1.0001, 9) {
//!     spline1.eval(x)?;
//!     spline2.eval(x)?;
//! }
//!
//! println!("{:?}", spline1.acc);
//! println!("{:?}", spline2.acc);
//! # Ok(())
//! # }
//! ```
//!
//! [`rust-GSL`]: https://lib.rs/crates/gsl
//! [`Accelerator`]: https://www.gnu.org/software/gsl/doc/html/interp.html#d-index-look-up-and-acceleration

mod acc;
mod error;
mod evals;
mod interp_types;
mod spline;

// Types refering to rgsl's Interp and Spline types, to avoid confusion.
pub(crate) type RgslSpline = rgsl::Spline;
pub(crate) type RgslInterpType = rgsl::InterpType;
pub(crate) type RgslInterpAccel = rgsl::InterpAccel;

pub use acc::Accelerator;
pub use error::SplineError;
pub use interp_types::InterpolationType;
pub use spline::Spline;

pub type Result<T, E = SplineError> = std::result::Result<T, E>;

#[allow(deprecated)]
pub(crate) type RgslValue = rgsl::Value;
