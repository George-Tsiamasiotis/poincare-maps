use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline};
use utils::array1D_getter_impl;

use crate::Result;

use ndarray::Array1;

/// q-factor reconstructed from a netCDF file.
pub struct Qfactor {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// Interpolation type.
    pub typ: String,

    /// Spline over the q-factor data, as a function of ψ_p.
    pub q_spline: DynSpline<f64>,
    /// Spline over the toroidal flux data, as a function of ψ_p.
    pub psi_spline: DynSpline<f64>,
    /// The value of the poloidal angle ψ_p at the wall.

    /// The value of the poloidal angle ψp at the wall.
    pub psip_wall: f64,
    /// The value of the toroidal angle ψ at the wall.
    pub psi_wall: f64,
}

/// Creation
impl Qfactor {
    /// Constructs a [`Qfactor`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Note
    ///
    /// The value `ψ = 0.0` is prepended at the ψ data array, and the first value of the q array is
    /// prepended (duplicated) in the q array, to assure correct interpolation near the magnetic axis.
    ///
    /// # Example
    /// ```no_run
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

        // Make path absolute. Just unwrap, Equilibrium checks if it exists.
        let path = std::path::absolute(path).unwrap();

        let eq = Equilibrium::from_file(&path)?;

        let psip_data = extract_1d_var(&eq.file, PSIP_COORD)?
            .as_standard_layout()
            .to_vec();
        let psi_data = extract_1d_var(&eq.file, PSI_COORD)?
            .as_standard_layout()
            .to_vec();
        let q_data = extract_1d_var(&eq.file, Q_FACTOR)?
            .as_standard_layout()
            .to_vec();

        let q_spline = make_spline(typ, &psip_data, &q_data)?;
        let psi_spline = make_spline(typ, &psip_data, &psi_data)?;

        let psip_wall = psip_data.last().copied().unwrap();
        let psi_wall = psi_data.last().copied().unwrap();

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            q_spline,
            psi_spline,
            psip_wall,
            psi_wall,
        })
    }
}

/// Interpolation
impl Qfactor {
    /// Calculates the q-factor `q(ψ_p)`.
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
    pub fn q(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.q_spline.eval(psip, acc)?)
    }

    /// Calculates the toroidal flux `ψ(ψ_p)`.
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
    pub fn psi(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.psi_spline.eval(psip, acc)?)
    }
}

impl Qfactor {
    /// Returns the `q` data calculated from `dψ/dψp` as a 1D array.
    pub fn q_data_derived(&self) -> Array1<f64> {
        let mut acc = Accelerator::new();
        self.q_spline
            .xa
            .iter()
            .map(|psip| self.psi_spline.eval_deriv(*psip, &mut acc).unwrap())
            .collect::<Vec<f64>>()
            .into()
    }
}

array1D_getter_impl!(Qfactor, psip_data, q_spline.xa);
array1D_getter_impl!(Qfactor, q_data, q_spline.ya);
array1D_getter_impl!(Qfactor, psi_data, psi_spline.ya);

impl std::fmt::Debug for Qfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Qfactor")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall))
            .field("ψ_wall", &format!("{:.7}", self.psi_wall))
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

        assert_eq!(q.psip_data().shape(), [101]);
        assert_eq!(q.q_data().shape(), [101]);
        assert_eq!(q.psi_data().shape(), [101]);
        assert_eq!(q.q_data_derived().shape(), [101]);
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
