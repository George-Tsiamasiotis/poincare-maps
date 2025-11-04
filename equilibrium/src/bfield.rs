use std::path::PathBuf;

use ndarray::{Array1, Array2};
use rsl_interpolation::{Accelerator, Cache, DynSpline2d};
use utils::array1D_getter_impl;

use crate::Result;

/// Magnetic field reconstructed from a netCDF file.
pub struct Bfield {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// Interpolation type.
    pub typ: String,

    /// Spline over the magnetic field strength data, as a function of ψp, θ.
    pub b_spline: DynSpline2d<f64>,
    /// Spline over the R coordinate, as a function of ψp, θ.
    pub r_spline: DynSpline2d<f64>,
    /// Spline over the Z coordinate, as a function of ψp, θ.
    pub z_spline: DynSpline2d<f64>,

    /// Magnetic field strength on the axis **in \[T\]**.
    pub baxis: f64,
    /// The tokamak's major radius **in \[m\]**.
    pub raxis: f64,
    /// The value of the poloidal angle ψp at the wall.
    pub psip_wall: f64,
    /// The value of the toroidal angle ψ at the wall.
    pub psi_wall: f64,
}

/// Creation
impl Bfield {
    /// Constructs a [`Bfield`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
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
        // For `psi_wall`
        let psi_data = extract_1d_var(&eq.file, PSI_COORD)?
            .as_standard_layout()
            .to_vec();
        let theta_data = eq.get_1d(THETA_COORD)?.to_vec();

        let b_data = eq.get_2d(B_FIELD)?;
        let r_data = eq.get_2d(R)?;
        let z_data = eq.get_2d(Z)?;
        let baxis = eq.get_scalar(B_AXIS)?;
        let raxis = eq.get_scalar(R_AXIS)?;
        // let zaxis_val = eq.get_scalar(Z_AXIS)?;

        // Transpose of gcmotion
        // let b_axis_values = Array2::from_elem((1, b_data.ncols()), 1.0); // B0 = 1 [NU]
        // let b_data = concatenate![Axis(0), b_axis_values, b_data]; // e.g. [101, 3620]

        // `Spline.za` is in Fortran order.
        let order = ndarray::Order::ColumnMajor;
        let b_data_flat = b_data.flatten_with_order(order).to_vec();

        // let r_axis_values = Array2::from_elem((1, r_data.ncols()), raxis_val);
        // let r_data = concatenate![Axis(0), r_axis_values, r_data];
        let r_data_flat = r_data.flatten_with_order(order).to_vec();
        //
        // let z_axis_values = Array2::from_elem((1, z_data.ncols()), zaxis_val);
        // let z_data = concatenate![Axis(0), z_axis_values, z_data];
        let z_data_flat = z_data.flatten_with_order(order).to_vec();

        let b_spline = make_spline2d(typ, &psip_data, &theta_data, &b_data_flat)?;
        let r_spline = make_spline2d(typ, &psip_data, &theta_data, &r_data_flat)?;
        let z_spline = make_spline2d(typ, &psip_data, &theta_data, &z_data_flat)?;

        let psip_wall = psip_data.last().copied().unwrap();
        let psi_wall = psi_data.last().copied().unwrap();

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            b_spline,
            r_spline,
            z_spline,
            baxis,
            raxis,
            psip_wall,
            psi_wall,
        })
    }
}

/// Interpolation
impl Bfield {
    /// Calculates `B(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let b =  bfield.b(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn b(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self.b_spline.eval(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕B(ψp, θ) /𝜕𝜃`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let db_dtheta = bfield.db_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn db_dtheta(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        // Ok(self.db_dtheta_spline.eval(psi, theta, xacc, yacc)?)
        Ok(self
            .b_spline
            .eval_deriv_y(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕B(ψp, θ) /𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let db_dpsip = bfield.db_dpsip(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn db_dpsip(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_x(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕²B(ψp, θ) /𝜕𝜓p²`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dpsip2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dpsip2(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_xx(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕²B(ψp, θ) /𝜕θ²`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dtheta2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dtheta2(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_yy(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕²B(ψp, θ) /𝜕ψp𝜕θ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dpsip_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dpsip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_xy(psip, mod2pi(theta), xacc, yacc, cache)?)
    }
}

/// Generates getters that return [T] fields to Array1<T>.
///
/// Since the data is internally stored as a Vec in Fortran order, we need a specialized macro.
///
/// This is needed for implementing python getter wrappers.
macro_rules! array2D_getter_impl {
    ($fun_name:ident, $spline:ident) => {
        pub fn $fun_name(&self) -> Array2<f64> {
            // `Spline.za` is in Fortran order.
            let shape = (self.$spline.ya.len(), self.$spline.xa.len());
            Array2::from_shape_vec(shape, self.$spline.za.to_vec())
                .expect("Error extracting 2D spline data")
                .reversed_axes()
        }
    };
}

/// Data Extraction
impl Bfield {
    array2D_getter_impl!(b_data, b_spline);
    array2D_getter_impl!(r_data, b_spline);
    array2D_getter_impl!(z_data, b_spline);

    /// Returns the `𝜕B(ψp, θ) /𝜕ψp` data as a 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    pub fn db_dpsip_data(&self) -> Array2<f64> {
        // `Spline.za` is in Fortran order.
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();

        let shape = (self.b_spline.ya.len(), self.b_spline.xa.len());

        let mut db_dpsip_vec = Vec::<f64>::with_capacity(shape.0 * shape.1);
        for j in 0..shape.0 {
            for i in 0..shape.1 {
                let psip = self.b_spline.xa[i];
                let theta = self.b_spline.ya[j];
                let db_dtheta = self
                    .db_dpsip(psip, theta, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                db_dpsip_vec.push(db_dtheta);
            }
        }

        let db_dpsip_grid = Array2::from_shape_vec(shape, db_dpsip_vec).unwrap();
        db_dpsip_grid.reversed_axes()
    }

    /// Returns the `𝜕B(ψp, θ) /𝜕𝜃` data as a 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    pub fn db_dtheta_data(&self) -> Array2<f64> {
        // `Spline.za` is in Fortran order.
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        let shape = (self.b_spline.ya.len(), self.b_spline.xa.len());

        let mut db_dtheta_vec = Vec::<f64>::with_capacity(shape.0 * shape.1);
        for j in 0..shape.0 {
            for i in 0..shape.1 {
                let psip = self.b_spline.xa[i];
                let theta = self.b_spline.ya[j];
                let db_dtheta = self
                    .db_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                db_dtheta_vec.push(db_dtheta);
            }
        }

        let db_dtheta_grid = Array2::from_shape_vec(shape, db_dtheta_vec).unwrap();
        db_dtheta_grid.reversed_axes()
    }
}

array1D_getter_impl!(Bfield, psip_data, b_spline.xa);
array1D_getter_impl!(Bfield, theta_data, b_spline.ya);

/// Returns θ % 2π.
fn mod2pi(theta: f64) -> f64 {
    use std::f64::consts::TAU;
    theta.rem_euclid(TAU)
}

impl std::fmt::Debug for Bfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Bfield")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("Baxis [T]", &format!("{:.7}", self.baxis))
            .field("Raxis [m]", &format!("{:.7}", self.raxis))
            .field("ψp_wall", &format!("{:.7}", self.psip_wall))
            .field("ψ_wall", &format!("{:.7}", self.psi_wall))
            .field("shape", &(self.b_spline.xa.len(), self.b_spline.ya.len()))
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn create_bfield() -> Bfield {
        let path = PathBuf::from("../data.nc");
        Bfield::from_dataset(&path, "bicubic").unwrap()
    }

    #[test]
    fn test_bfield_creation() {
        create_bfield();
    }

    #[test]
    fn test_extraction_methods() {
        let b = create_bfield();
        let _ = format!("{b:?}");

        assert_eq!(b.psip_data().ndim(), 1);
        assert_eq!(b.theta_data().ndim(), 1);
        assert_eq!(b.r_data().ndim(), 2);
        assert_eq!(b.z_data().ndim(), 2);
        assert_eq!(b.b_data().ndim(), 2);
        assert_eq!(b.db_dpsip_data().ndim(), 2);
        assert_eq!(b.db_dtheta_data().ndim(), 2);
    }

    #[test]
    fn test_spline_evaluation() {
        let b = create_bfield();
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();

        let psip = 0.015;
        let theta = 0.0;
        b.b(psip, theta, &mut xacc, &mut yacc, &mut cache).unwrap();
        b.db_dpsip(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.db_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.d2b_dpsip2(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.d2b_dtheta2(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.d2b_dpsip_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
    }
}
