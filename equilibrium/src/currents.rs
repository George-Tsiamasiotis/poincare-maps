use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline, make_spline};
use utils::array1D_getter_impl;

use crate::Flux;
use crate::Result;

use ndarray::Array1;
use safe_unwrap::safe_unwrap;

/// Plasma current reconstructed from a netCDF file.
pub struct Currents {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    pub typ: String,
    /// Spline over the g-current data, as a function of ψp.
    pub g_spline: DynSpline<f64>,
    /// Spline over the I-current data, as a function of ψp.
    pub i_spline: DynSpline<f64>,
}

/// Creation
impl Currents {
    /// Constructs a [`Currents`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str) -> Result<Self> {
        use crate::extract::*;
        use config::netcdf_fields::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP)?.as_standard_layout().to_owned();
        let g_data = extract_1d_array(&f, G)?.as_standard_layout().to_owned();
        let i_data = extract_1d_array(&f, I)?.as_standard_layout().to_owned();

        let g_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", g_data.as_slice()),
        )?;
        let i_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", i_data.as_slice()),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            g_spline,
            i_spline,
        })
    }
}

// Interpolation
impl Currents {
    /// Calculates `g(ψp)`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let g = currents.g(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn g(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.g_spline.eval(psip, acc)?)
    }

    /// Calculates `I(ψp)`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let i = currents.i(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn i(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.i_spline.eval(psip, acc)?)
    }

    /// Calculates `𝜕g(ψp)/𝜕ψp`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dg = currents.dg_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dg_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.g_spline.eval_deriv(psip, acc)?)
    }

    /// Calculates `𝜕I(ψp)/𝜕ψp`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let di = currents.di_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn di_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.i_spline.eval_deriv(psip, acc)?)
    }
}

// Data extraction
impl Currents {
    array1D_getter_impl!(psip_data, g_spline.xa, Flux);
    array1D_getter_impl!(g_data, g_spline.ya, f64);
    array1D_getter_impl!(i_data, i_spline.ya, f64);

    /// Returns the value of the poloidal angle ψp at the wall.
    pub fn psip_wall(&self) -> Flux {
        safe_unwrap!("ya is non-empty", self.g_spline.xa.last().copied())
    }
}

impl std::fmt::Debug for Currents {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Current")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall()))
            .field("len", &self.g_spline.xa.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use config::STUB_NETCDF_PATH;

    fn create_current() -> Currents {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        Currents::from_dataset(&path, "akima").unwrap()
    }

    #[test]
    fn test_current_creation() {
        create_current();
    }

    #[test]
    fn test_extraction_methods() {
        let c = create_current();
        let _ = format!("{c:?}");

        assert_eq!(c.psip_data().ndim(), 1);
        assert_eq!(c.g_data().ndim(), 1);
        assert_eq!(c.i_data().ndim(), 1);
    }

    #[test]
    fn test_spline_evaluation() {
        let c = create_current();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        c.g(psip, &mut acc).unwrap();
        c.i(psip, &mut acc).unwrap();
        c.di_dpsip(psip, &mut acc).unwrap();
        c.dg_dpsip(psip, &mut acc).unwrap();
    }
}
