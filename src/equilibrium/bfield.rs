use std::path::PathBuf;

use ndarray::concatenate;
use ndarray::{Array1, Array2, Axis};
use numpy::{PyArray1, PyArray2, ToPyArray};
use rsl_interpolation::{Accelerator, Cache, DynSpline2d};

use crate::Result;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Magnetic field reconstructed from a netCDF file.
#[pyclass(frozen, immutable_type)]
pub struct Bfield {
    /// Path to the netCDF file.
    #[pyo3(get)]
    pub path: PathBuf,
    /// Interpolation type.
    #[pyo3(get)]
    pub typ: String,

    /// Spline over the magnetic field strength data, as a function of ψp, θ.
    pub b_spline: DynSpline2d<f64>,
    /// Spline over the R coordinate, as a function of ψp, θ.
    pub r_spline: DynSpline2d<f64>,
    /// Spline over the Z coordinate, as a function of ψp, θ.
    pub z_spline: DynSpline2d<f64>,

    /// Magnetic field strength on the axis **in \[T\]**.
    #[pyo3(get)]
    pub baxis: f64,
    /// The tokamak's major radius **in \[m\]**.
    #[pyo3(get)]
    pub raxis: f64,
    /// The value of the poloidal angle ψp at the wall.
    #[pyo3(get)]
    pub psip_wall: f64,
    /// The value of the toroidal angle ψ at the wall.
    #[pyo3(get)]
    pub psi_wall: f64,
}

/// Wrapper methods exposed to Python.
#[pymethods]
impl Bfield {
    /// Creates a new [`Bfield`]
    ///
    /// Wrapper around [`Bfield::from_dataset`]. This is a workaround to return a [`PyErr`].
    #[coverage(off)]
    #[new]
    pub fn new_py(path: &str, typ: &str) -> PyResult<Self> {
        let path = PathBuf::from(path);
        match Self::from_dataset(&path, typ) {
            Ok(bfield) => Ok(bfield),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    /// Returns the `psip` coordinate data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "psip_data")]
    pub fn psip_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.psip_data().to_pyarray(py)
    }

    /// Returns the `theta` coordinate data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "theta_data")]
    pub fn theta_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.theta_data().to_pyarray(py)
    }

    /// Returns the `B` data grid as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "b_data")]
    pub fn b_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.b_data().to_pyarray(py)
    }

    /// Returns the `R` data grid as  Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "r_data")]
    pub fn r_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.r_data().to_pyarray(py)
    }

    /// Returns the `Z` data grid as  Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "z_data")]
    pub fn z_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.z_data().to_pyarray(py)
    }

    /// Returns the `𝜕B(ψp, θ) /𝜕ψp` data as a Numpy 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    #[pyo3(name = "db_dpsip_data")]
    #[coverage(off)]
    pub fn db_dspip_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.db_dpsip_data().to_pyarray(py)
    }

    /// Returns the `𝜕B(ψp, θ) /𝜕𝜃` data as a Numpy 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    #[pyo3(name = "db_dtheta_data")]
    #[coverage(off)]
    pub fn db_dtheta_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.db_dtheta_data().to_pyarray(py)
    }

    #[coverage(off)]
    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }
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

        // Make path absolute. Just unwrap, Equilibrium checks if it exists.
        let path = std::path::absolute(path).unwrap();

        let eq = Equilibrium::from_file(&path)?;

        // Add 0.0 manualy, which corresponds to the axis value.
        let psip_data = extract_var_with_axis_value(&eq.file, PSIP_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        // For `psi_wall`
        let psi_data = extract_var_with_axis_value(&eq.file, PSI_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        let theta_data = eq.get_1d(THETA_COORD)?.to_vec();

        let b_data = eq.get_2d(B_FIELD)?;
        let r_data = eq.get_2d(R)?;
        let z_data = eq.get_2d(Z)?;
        let baxis_val = eq.get_scalar(B_AXIS)?;
        let raxis_val = eq.get_scalar(R_AXIS)?;
        let zaxis_val = eq.get_scalar(Z_AXIS)?;

        // Transpose of gcmotion
        let b_axis_values = Array2::from_elem((1, b_data.ncols()), 1.0); // B0 = 1 [NU]
        let b_data = concatenate![Axis(0), b_axis_values, b_data]; // e.g. [101, 3620]

        // `Spline.za` is in Fortran order.
        let order = ndarray::Order::ColumnMajor;
        let b_data_flat = b_data.flatten_with_order(order).to_vec();

        let r_axis_values = Array2::from_elem((1, r_data.ncols()), raxis_val);
        let r_data = concatenate![Axis(0), r_axis_values, r_data];
        let r_data_flat = r_data.flatten_with_order(order).to_vec();

        let z_axis_values = Array2::from_elem((1, z_data.ncols()), zaxis_val);
        let z_data = concatenate![Axis(0), z_axis_values, z_data];
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
            baxis: baxis_val,
            raxis: raxis_val,
            psip_wall,
            psi_wall,
        })
    }

    /// Returns the `psip` coordinate data as a 1D array.
    pub fn psip_data(&self) -> Array1<f64> {
        Array1::from_vec(self.b_spline.xa.to_vec())
    }

    /// Returns the `theta` coordinate data as a 1D array.
    pub fn theta_data(&self) -> Array1<f64> {
        Array1::from_vec(self.b_spline.ya.to_vec())
    }

    /// Returns the `R` data grid as a 2D array.
    pub fn r_data(&self) -> Array2<f64> {
        let shape = (self.r_spline.ya.len(), self.r_spline.xa.len());
        match Array2::from_shape_vec(shape, self.r_spline.za.to_vec()) {
            Ok(r_grid) => r_grid.reversed_axes(),
            Err(_) => unreachable!(),
        }
    }

    /// Returns the `Z` data grid as a 2D array.
    pub fn z_data(&self) -> Array2<f64> {
        let shape = (self.z_spline.ya.len(), self.z_spline.xa.len());
        match Array2::from_shape_vec(shape, self.z_spline.za.to_vec()) {
            Ok(z_grid) => z_grid.reversed_axes(),
            Err(_) => unreachable!(),
        }
    }

    /// Returns the `B` data grid as a 2D array.
    pub fn b_data(&self) -> Array2<f64> {
        let shape = (self.b_spline.ya.len(), self.b_spline.xa.len());
        match Array2::from_shape_vec(shape, self.b_spline.za.to_vec()) {
            // `Spline.za` is in Fortran order.
            Ok(b_grid) => b_grid.reversed_axes(),
            Err(_) => unreachable!(),
        }
    }

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

impl Bfield {
    /// Calculates `B(ψp, θ)`,
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
        let path = PathBuf::from("./data.nc");
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

        assert_eq!(b.psip_data().shape(), [101]);
        assert_eq!(b.theta_data().shape(), [3620]);
        assert_eq!(b.r_data().shape(), [101, 3620]);
        assert_eq!(b.z_data().shape(), [101, 3620]);
        assert_eq!(b.b_data().shape(), [101, 3620]);
        assert_eq!(b.db_dpsip_data().shape(), [101, 3620]);
        assert_eq!(b.db_dtheta_data().shape(), [101, 3620]);
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
