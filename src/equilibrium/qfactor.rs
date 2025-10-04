use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline};

use crate::Result;

/// q-factor reconstructed from a netCDF file.
#[non_exhaustive]
pub struct Qfactor {
    /// Spline over the q-factor data, as a function of ψ_p.
    pub q_spline: DynSpline<f64>,
    /// Spline over the toroidal flux data, as a function of ψ_p.
    pub psi_spline: DynSpline<f64>,
}

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
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        let eq = Equilibrium::from_file(path)?;

        // Add 0.0 manualy, which corresponds to q0.
        let psip_data = extract_var_with_axis_value(&eq.file, PSIP_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        let psi_data = extract_var_with_axis_value(&eq.file, PSI_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        // Manually add q0 to the array.
        let q_data = extract_var_with_first_axis_value(&eq.file, Q_FACTOR)?
            .as_standard_layout()
            .to_vec();

        let q_spline = make_spline(typ, &psip_data, &q_data)?;
        let psi_spline = make_spline(typ, &psip_data, &psi_data)?;

        Ok(Self {
            q_spline,
            psi_spline,
        })
    }
}

impl Qfactor {
    /// Calculates the q-factor `q(ψ_p)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let q =  qfactor.q(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn q(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(psip.is_sign_positive());
        Ok(self.q_spline.eval(psip, acc)?)
    }

    /// Calculates the toroidal flux `ψ(ψ_p)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let psi =  qfactor.psi(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn psi(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(psip.is_sign_positive());
        Ok(self.psi_spline.eval(psip, acc)?)
    }
}
