use std::path::PathBuf;

use ndarray::concatenate;
use ndarray::{Array2, Axis};
use rsl_interpolation::{Accelerator, DynSpline2d};

use crate::Result;

/// Magnetic field reconstructed from a netCDF file.
pub struct Bfield {
    /// Spline over the magnetic field strength data, as a function of ψ_p.
    b_spline: DynSpline2d<f64>,
}

impl Bfield {
    /// Constructs a [`Bfield`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        let eq = Equilibrium::from_file(path)?;

        // Add 0.0 manualy, which corresponds to the axis value.
        let psip_data = extract_var_with_axis_value(&eq.file, PSIP_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        let theta_data = eq.get_1d(THETA_COORD)?.to_vec();

        let b_data = eq.get_2d(B_FIELD)?;

        // Transpose of gcmotion
        let b_axis_values = Array2::from_elem((1, b_data.ncols()), 1.0); // B0 = 1 [NU]
        let b_data = concatenate![Axis(0), b_axis_values, b_data]; // e.g. [101, 3620]
        let b_data_flat = b_data.flatten().to_vec();

        let b_spline = make_spline2d(typ, &psip_data, &theta_data, &b_data_flat)?;

        Ok(Self { b_spline })
    }
}

impl Bfield {
    /// Calculates `B(ψ_p, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let b =  bfield.b(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn b(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<f64> {
        Ok(self.b_spline.eval(psip, mod_theta(theta), xacc, yacc)?)
    }

    /// Calculates `𝜕B(ψ_p, θ) /𝜕𝜃`.
    ///
    /// # Example
    ///
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let db_dtheta = bfield.db_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn db_dtheta(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<f64> {
        // Ok(self.db_dtheta_spline.eval(psi, theta, xacc, yacc)?)
        Ok(self
            .b_spline
            .eval_deriv_y(psip, mod_theta(theta), xacc, yacc)?)
    }

    /// Calculates `𝜕B(ψ_p, θ) /𝜕ψ_p`.
    ///
    /// # Example
    ///
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let db_dpsi = bfield.db_dpsi(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn db_dpsi(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<f64> {
        // Ok(self.db_dpsi_spline.eval(psi, theta, xacc, yacc)?)
        Ok(self
            .b_spline
            .eval_deriv_x(psip, mod_theta(theta), xacc, yacc)?)
    }

    /// Calculates `𝜕²B(ψ_p, θ) /𝜕𝜓_p²`.
    ///
    /// # Example
    ///
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let d2b_dpsi2 = bfield.d2b_dpsi2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dpsi2(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_xx(psip, mod_theta(theta), xacc, yacc)?)
    }
}

/// Returns θ % 2π.
fn mod_theta(theta: f64) -> f64 {
    use std::f64::consts::TAU;
    theta.rem_euclid(TAU)
}
