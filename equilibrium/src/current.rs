use std::path::PathBuf;

use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynSpline};
use utils::array1D_getter_impl;

use crate::Result;

/// Plasma current reconstructed from a netCDF file.
pub struct Current {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// Interpolation type.
    pub typ: String,

    /// Spline over the g-current data, as a function of ψp.
    pub g_spline: DynSpline<f64>,
    /// Spline over the I-current data, as a function of ψp.
    pub i_spline: DynSpline<f64>,

    /// The value of the poloidal angle ψp at the wall.
    pub psip_wall: f64,
}

/// Creation
impl Current {
    /// Constructs a [`Current`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Note
    ///
    /// The value `ψ = 0.0` is prepended at the ψ data array, and the first values of the i and g arrays
    /// is prepended (duplicated) in each array, to assure correct interpolation near the magnetic axis.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let cur = Current::from_dataset(&path, "cubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        // Make path absolute. Just unwrap, Equilibrium checks if it exists.
        let path = std::path::absolute(path).unwrap();

        let eq = Equilibrium::from_file(&path)?;

        let psip_data = extract_1d_var(&eq.file, PSIP_COORD)?
            .as_standard_layout()
            .to_vec();
        let g_data = extract_1d_var(&eq.file, CURRENT_G)?
            .as_standard_layout()
            .to_vec();
        let i_data = extract_1d_var(&eq.file, CURRENT_I)?
            .as_standard_layout()
            .to_vec();

        let g_spline = make_spline(typ, &psip_data, &g_data)?;
        let i_spline = make_spline(typ, &psip_data, &i_data)?;

        let psip_wall = psip_data.last().copied().unwrap();

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            g_spline,
            i_spline,
            psip_wall,
        })
    }
}

/// Interpolation
impl Current {
    /// Calculates `g(ψp)`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let current = Current::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let g = current.g(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn g(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
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
    /// let path = PathBuf::from("../data.nc");
    /// let current = Current::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let i = current.i(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn i(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
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
    /// let path = PathBuf::from("../data.nc");
    /// let current = Current::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dg = current.dg_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
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
    /// let path = PathBuf::from("../data.nc");
    /// let current = Current::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let di = current.di_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.i_spline.eval_deriv(psip, acc)?)
    }
}

array1D_getter_impl!(Current, psip_data, g_spline.xa);
array1D_getter_impl!(Current, g_data, g_spline.ya);
array1D_getter_impl!(Current, i_data, g_spline.ya);

impl std::fmt::Debug for Current {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Current")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall))
            .field("len", &self.g_spline.xa.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn create_current() -> Current {
        let path = PathBuf::from("../data.nc");
        Current::from_dataset(&path, "akima").unwrap()
    }

    #[test]
    fn test_current_creation() {
        create_current();
    }

    #[test]
    fn test_extraction_methods() {
        let c = create_current();
        let _ = format!("{c:?}");

        assert_eq!(c.psip_data().shape(), [101]);
        assert_eq!(c.g_data().shape(), [101]);
        assert_eq!(c.i_data().shape(), [101]);
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
