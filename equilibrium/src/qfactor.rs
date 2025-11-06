use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline};
use utils::array1D_getter_impl;

use crate::Flux;
use crate::Result;

use ndarray::Array1;
use safe_unwrap::safe_unwrap;

/// q-factor reconstructed from a netCDF file.
pub struct Qfactor {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// Interpolation type.
    pub typ: String,
    /// Spline over the q-factor data, as a function of ψp.
    pub q_spline: DynSpline<f64>,
    /// Spline over the toroidal flux data, as a function of ψp.
    pub psi_spline: DynSpline<f64>,
}

// Creation
impl Qfactor {
    /// Constructs a [`Qfactor`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;

        let eq = Equilibrium::from_file(&path)?;

        let psip_data = extract_1d_var(&eq.file, PSIP_COORD)?
            .as_standard_layout()
            .to_owned();
        let psi_data = extract_1d_var(&eq.file, PSI_COORD)?
            .as_standard_layout()
            .to_owned();
        let q_data = extract_1d_var(&eq.file, Q_FACTOR)?
            .as_standard_layout()
            .to_owned();

        let q_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", q_data.as_slice()),
        )?;
        let psi_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", psi_data.as_slice()),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            q_spline,
            psi_spline,
        })
    }
}

// Interpolation
impl Qfactor {
    /// Calculates the q-factor `q(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let q =  qfactor.q(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn q(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.q_spline.eval(psip, acc)?)
    }

    /// Calculates the toroidal flux `ψ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let psi =  qfactor.psi(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn psi(&self, psip: Flux, acc: &mut Accelerator) -> Result<Flux> {
        Ok(self.psi_spline.eval(psip, acc)?)
    }
}

// Data extraction
impl Qfactor {
    array1D_getter_impl!(psip_data, q_spline.xa, Flux);
    array1D_getter_impl!(psi_data, psi_spline.ya, Flux);
    array1D_getter_impl!(q_data, q_spline.ya, f64);

    /// Returns the `q` data calculated from `dψ/dψp` as a 1D array.
    pub fn q_data_derived(&self) -> Result<Array1<f64>> {
        let mut acc = Accelerator::new();
        let mut q_data = Array1::from_elem(self.q_spline.xa.len(), f64::NAN);

        for (i, psip) in self.psip_data().iter().enumerate() {
            q_data[[i]] = self.psi_spline.eval_deriv(*psip, &mut acc)?;
        }

        Ok(q_data)
    }

    /// Returns the value of the poloidal angle ψp at the wall.
    pub fn psip_wall(&self) -> Flux {
        safe_unwrap!("ya is non-empty", self.q_spline.xa.last().copied())
    }

    /// Returns the value of the toroidal angle ψ at the wall.
    pub fn psi_wall(&self) -> Flux {
        safe_unwrap!("xa is non-empty", self.psi_spline.ya.last().copied())
    }
}

impl std::fmt::Debug for Qfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Qfactor")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall()))
            .field("ψ_wall", &format!("{:.7}", self.psi_wall()))
            .field("len", &self.q_spline.xa.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn create_qfactor() -> Qfactor {
        let path = PathBuf::from("../data.nc");
        Qfactor::from_dataset(&path, "akima").unwrap()
    }

    #[test]
    fn test_qfactor_creation() {
        create_qfactor();
    }

    #[test]
    fn test_extraction_methods() {
        let q = create_qfactor();
        let _ = format!("{q:?}");

        assert_eq!(q.psip_data().ndim(), 1);
        assert_eq!(q.q_data().ndim(), 1);
        assert_eq!(q.psi_data().ndim(), 1);
        assert_eq!(q.q_data_derived().unwrap().ndim(), 1);
    }

    #[test]
    fn test_spline_evaluation() {
        let q = create_qfactor();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        q.q(psip, &mut acc).unwrap();
        q.psi(psip, &mut acc).unwrap();
    }
}
