use std::path::PathBuf;

use rsl_interpolation::{Accelerator, Cache, DynSpline2d, make_spline2d};
use utils::array1D_getter_impl;

use crate::Result;
use crate::{Flux, Radians};

use ndarray::{Array1, Array2};
use safe_unwrap::safe_unwrap;

/// Magnetic field reconstructed from a netCDF file.
pub struct Bfield {
    /// Path to the netCDF file.
    pub path: PathBuf,
    /// 2D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.Interp2dType.html#implementors
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
}

// Creation
impl Bfield {
    /// Constructs a [`Bfield`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
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
        let theta_data = extract_1d_array(&f, THETA)?.as_standard_layout().to_owned();

        let b_data = extract_2d_array(&f, B)?;
        let r_data = extract_2d_array(&f, R)?;
        let z_data = extract_2d_array(&f, Z)?;
        let baxis = extract_scalar(&f, BAXIS)?;
        let raxis = extract_scalar(&f, RAXIS)?;

        // `Spline.za` is in Fortran order.
        let order = ndarray::Order::ColumnMajor;
        let b_data_flat = b_data.flatten_with_order(order).to_owned();
        let r_data_flat = r_data.flatten_with_order(order).to_owned();
        let z_data_flat = z_data.flatten_with_order(order).to_owned();

        let b_spline = make_spline2d(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", theta_data.as_slice()),
            safe_unwrap!("array is non-empty", b_data_flat.as_slice()),
        )?;
        let r_spline = make_spline2d(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", theta_data.as_slice()),
            safe_unwrap!("array is non-empty", r_data_flat.as_slice()),
        )?;
        let z_spline = make_spline2d(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", theta_data.as_slice()),
            safe_unwrap!("array is non-empty", z_data_flat.as_slice()),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            b_spline,
            r_spline,
            z_spline,
            baxis,
            raxis,
        })
    }
}

// Interpolation
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
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
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self.b_spline.eval(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕B(ψp, θ) /𝜕𝜃`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
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
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_y(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕B(ψp, θ) /𝜕ψp`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
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
        psip: Flux,
        theta: Radians,
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
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
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
        psip: Flux,
        theta: Radians,
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
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
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
        psip: Flux,
        theta: Radians,
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
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
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
        psip: Flux,
        theta: Radians,
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
            safe_unwrap!(
                "checked by spline",
                Array2::from_shape_vec(shape, self.$spline.za.to_vec())
            )
            .reversed_axes()
        }
    };
}

// Data Extraction
impl Bfield {
    array1D_getter_impl!(psip_data, b_spline.xa, Flux);
    array1D_getter_impl!(theta_data, b_spline.ya, Flux);
    array2D_getter_impl!(b_data, b_spline);
    array2D_getter_impl!(r_data, r_spline);
    array2D_getter_impl!(z_data, z_spline);

    /// Returns the `𝜕B(ψp, θ) /𝜕ψp` data as a 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    pub fn db_dpsip_data(&self) -> Result<Array2<f64>> {
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
                let db_dtheta = self.db_dpsip(psip, theta, &mut xacc, &mut yacc, &mut cache)?;
                db_dpsip_vec.push(db_dtheta);
            }
        }

        let db_dpsip_grid = Array2::from_shape_vec(shape, db_dpsip_vec)?;
        Ok(db_dpsip_grid.reversed_axes())
    }

    /// Returns the `𝜕B(ψp, θ) /𝜕𝜃` data as a 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    pub fn db_dtheta_data(&self) -> Result<Array2<f64>> {
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
                let db_dtheta = self.db_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)?;
                db_dtheta_vec.push(db_dtheta);
            }
        }

        let db_dtheta_grid = Array2::from_shape_vec(shape, db_dtheta_vec)?;
        Ok(db_dtheta_grid.reversed_axes())
    }

    /// Returns the value of the poloidal angle ψp at the wall.
    pub fn psip_wall(&self) -> Flux {
        safe_unwrap!("ya is non-empty", self.b_spline.xa.last().copied())
    }

    /// Returns the value of the magnetic field strength B0 at the axis **in \[T\]**.
    pub fn baxis(&self) -> f64 {
        safe_unwrap!("ya is non-empty", self.b_spline.ya.last().copied())
    }

    /// Returns the value of the major radius R **in \[m\]**.
    pub fn raxis(&self) -> f64 {
        self.raxis
    }
}

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
            .field("ψp_wall", &format!("{:.7}", self.psip_wall()))
            .field("shape", &(self.b_spline.xa.len(), self.b_spline.ya.len()))
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use config::STUB_NETCDF_PATH;

    fn create_bfield() -> Bfield {
        let path = PathBuf::from(STUB_NETCDF_PATH);
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
        let _: f64 = b.baxis();
        let _: f64 = b.raxis();

        assert_eq!(b.psip_data().ndim(), 1);
        assert_eq!(b.theta_data().ndim(), 1);
        assert_eq!(b.r_data().ndim(), 2);
        assert_eq!(b.z_data().ndim(), 2);
        assert_eq!(b.b_data().ndim(), 2);
        assert_eq!(b.db_dpsip_data().unwrap().ndim(), 2);
        assert_eq!(b.db_dtheta_data().unwrap().ndim(), 2);
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
