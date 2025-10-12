use std::path::PathBuf;

use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynSpline};

use crate::Result;

use numpy::{PyArray1, ToPyArray};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Plasma current reconstructed from a netCDF file.
#[pyclass(frozen, immutable_type)]
pub struct Current {
    /// Path to the netCDF file.
    #[pyo3(get)]
    pub path: PathBuf,
    /// Interpolation type.
    #[pyo3(get)]
    pub typ: String,

    /// Spline over the g-current data, as a function of ψp.
    pub g_spline: DynSpline<f64>,
    /// Spline over the I-current data, as a function of ψp.
    pub i_spline: DynSpline<f64>,

    /// The value of the poloidal angle ψp at the wall.
    #[pyo3(get)]
    pub psip_wall: f64,
    /// The value of the toroidal angle ψ at the wall.
    #[pyo3(get)]
    pub psi_wall: f64,
}

#[pymethods]
impl Current {
    /// Creates a new [`Current`]
    ///
    /// Wrapper around [`Current::from_dataset`]. This is a workaround to return a [`PyErr`].
    #[coverage(off)]
    #[new]
    pub fn new(path: &str, typ: &str) -> PyResult<Self> {
        let path = PathBuf::from(path);
        match Self::from_dataset(&path, typ) {
            Ok(current) => Ok(current),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    /// Returns the `psip` coordinate data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "psip_data")]
    pub fn psip_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.psip_data().to_pyarray(py)
    }

    /// Returns the `g` data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "g_data")]
    pub fn g_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.g_data().to_pyarray(py)
    }

    /// Returns the `I` data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "i_data")]
    pub fn i_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.i_data().to_pyarray(py)
    }

    /// Returns the `𝜕g(ψp)/𝜕ψp` data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "dg_dpsip_data")]
    pub fn dg_dpsip_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.dg_dpsip_data().to_pyarray(py)
    }

    /// Returns the `𝜕I(ψp)/𝜕ψp` data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "di_dpsip_data")]
    pub fn di_dpsip_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.di_dpsip_data().to_pyarray(py)
    }

    #[coverage(off)]
    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }
}

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
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
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

        // Add 0.0 manualy, which corresponds to q0.
        let psip_data = extract_var_with_axis_value(&eq.file, PSIP_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        // For `psi_wall`
        let psi_data = extract_var_with_axis_value(&eq.file, PSI_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        // Manually add q0 to the array.
        let g_data = extract_var_with_first_axis_value(&eq.file, CURRENT_G)?
            .as_standard_layout()
            .to_vec();
        let i_data = extract_var_with_first_axis_value(&eq.file, CURRENT_I)?
            .as_standard_layout()
            .to_vec();

        let g_spline = make_spline(typ, &psip_data, &g_data)?;
        let i_spline = make_spline(typ, &psip_data, &i_data)?;

        let psip_wall = psip_data.last().copied().unwrap();
        let psi_wall = psi_data.last().copied().unwrap();

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            g_spline,
            i_spline,
            psip_wall,
            psi_wall,
        })
    }

    /// Returns the `psip` coordinate data as a 1D array.
    pub fn psip_data(&self) -> Array1<f64> {
        Array1::from_vec(self.g_spline.xa.to_vec())
    }

    /// Returns the `g` data as a 1D array.
    pub fn g_data(&self) -> Array1<f64> {
        Array1::from_vec(self.g_spline.ya.to_vec())
    }

    /// Returns the `I` data as a 1D array.
    pub fn i_data(&self) -> Array1<f64> {
        Array1::from_vec(self.i_spline.ya.to_vec())
    }

    /// Returns the `𝜕g(ψp)/𝜕ψp` data as a 1D array.
    pub fn dg_dpsip_data(&self) -> Array1<f64> {
        let mut acc = Accelerator::new();
        let psip_data = &self.g_spline.xa;
        let dg_dpsip_vec: Vec<f64> = psip_data
            .iter()
            .map(|psip| self.dg_dpsip(*psip, &mut acc).unwrap())
            .collect();

        Array1::from_vec(dg_dpsip_vec)
    }

    /// Returns the `𝜕I(ψp)/𝜕ψp` data as a 1D array.
    pub fn di_dpsip_data(&self) -> Array1<f64> {
        let mut acc = Accelerator::new();
        let psip_data = &self.i_spline.xa;
        let di_dpsip_vec: Vec<f64> = psip_data
            .iter()
            .map(|psip| self.di_dpsip(*psip, &mut acc).unwrap())
            .collect();

        Array1::from_vec(di_dpsip_vec)
    }
}

impl Current {
    /// Calculates `g(ψp)`
    ///
    /// # Example
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
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
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
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
}

impl Current {
    /// Calculates `𝜕g(ψp)/𝜕ψp`
    ///
    /// # Example
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
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
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
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

impl std::fmt::Debug for Current {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Current")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall))
            .field("ψ_wall", &format!("{:.7}", self.psi_wall))
            .field("len", &self.g_spline.xa.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn create_current() -> Current {
        let path = PathBuf::from("./data.nc");
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
        assert_eq!(c.dg_dpsip_data().shape(), [101]);
        assert_eq!(c.di_dpsip_data().shape(), [101]);
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
