use std::path::PathBuf;

use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynSpline};

use crate::Result;

use numpy::{PyArray1, ToPyArray};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Plasma current reconstructed from a netCDF file.
#[pyclass]
pub struct Current {
    /// Path to the netCDF file.
    path: PathBuf,
    /// Interpolation type.
    typ: Box<str>,
    /// Spline over the g-current data, as a function of ψ_p.
    g_spline: DynSpline<f64>,
    /// Spline over the I-current data, as a function of ψ_p.
    i_spline: DynSpline<f64>,
}

#[pymethods]
impl Current {
    #[new]
    /// Wrapper around `Current::from_dataset`.
    ///
    /// This is a workaround to return a `PyErr`.
    pub fn new(path: &str, typ: &str) -> PyResult<Self> {
        let path = PathBuf::from(path);
        match Self::from_dataset(&path, typ) {
            Ok(current) => Ok(current),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    pub fn psip_data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.g_spline.xa.to_pyarray(py)
    }

    pub fn g_data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.g_spline.ya.to_pyarray(py)
    }

    pub fn i_data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.i_spline.ya.to_pyarray(py)
    }

    pub fn dg_dpsip_data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut acc = Accelerator::new();
        let psip_data = &self.g_spline.xa;
        let dg_dpsip_vec: Vec<f64> = psip_data
            .iter()
            .map(|psip| self.dg_dpsip(*psip, &mut acc).unwrap())
            .collect();

        let dg_dpsip = Array1::from_vec(dg_dpsip_vec);

        dg_dpsip.to_pyarray(py)
    }

    pub fn di_dpsip_data<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut acc = Accelerator::new();
        let psip_data = &self.g_spline.xa;
        let di_dpsip_vec: Vec<f64> = psip_data
            .iter()
            .map(|psip| self.di_dpsip(*psip, &mut acc).unwrap())
            .collect();

        let di_dpsip = Array1::from_vec(di_dpsip_vec);

        di_dpsip.to_pyarray(py)
    }

    fn __repr__(&self) -> String {
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

        let eq = Equilibrium::from_file(path)?;

        // Add 0.0 manualy, which corresponds to q0.
        let psip_data = extract_var_with_axis_value(&eq.file, PSIP_COORD, 0.0)?
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

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            g_spline,
            i_spline,
        })
    }
}

impl Current {
    /// Calculates `g(ψ_p)`
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

    /// Calculates `I(ψ_p)`
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
    /// Calculates `𝜕g(ψ_p)/𝜕ψ_p`
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

    /// Calculates `𝜕I(ψ_p)/𝜕ψ_p`
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
            .field("len", &self.g_spline.xa.len())
            .finish()
    }
}
