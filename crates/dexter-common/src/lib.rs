//! Common utilities across the workspace.

mod interp;
mod macros;
mod threads;

pub use interp::{
    DynInterpolator, DynInterpolator2d, Interpolation1dType, Interpolation2dType, make_interp,
    make_interp2d,
};
pub use threads::{get_max_threads, set_num_threads};
