use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline};

use crate::Result;

use ndarray::Array1;
use numpy::{PyArray1, ToPyArray};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// q-factor reconstructed from a netCDF file.
#[pyclass(frozen, immutable_type)]
pub struct Qfactor {
    /// Path to the netCDF file.
    #[pyo3(get)]
    pub path: PathBuf,
    /// Interpolation type.
    #[pyo3(get)]
    pub typ: String,

    /// Spline over the q-factor data, as a function of ψ_p.
    pub q_spline: DynSpline<f64>,
    /// Spline over the toroidal flux data, as a function of ψ_p.
    pub psi_spline: DynSpline<f64>,
    /// The value of the poloidal angle ψ_p at the wall.

    /// The value of the poloidal angle ψp at the wall.
    #[pyo3(get)]
    pub psip_wall: f64,
    /// The value of the toroidal angle ψ at the wall.
    #[pyo3(get)]
    pub psi_wall: f64,
}

#[pymethods]
impl Qfactor {
    /// Creates a new [`Qfactor`]
    ///
    /// Wrapper around [`Qfactor::from_dataset`]. This is a workaround to return a [`PyErr`].
    #[coverage(off)]
    #[new]
    pub fn new(path: &str, typ: &str) -> PyResult<Self> {
        let path = PathBuf::from(path);
        match Self::from_dataset(&path, typ) {
            Ok(qfactor) => Ok(qfactor),
            Err(err) => Err(PyTypeError::new_err(err.to_string())),
        }
    }

    /// Returns the `psip` coordinate data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "psip_data")]
    pub fn psip_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.psip_data().to_pyarray(py)
    }

    /// Returns the `q` data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "q_data")]
    pub fn q_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.q_data().to_pyarray(py)
    }

    /// Returns the `psi` data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "psi_data")]
    pub fn psi_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.psi_data().to_pyarray(py)
    }

    #[coverage(off)]
    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }
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

        // Make path absolute. Just unwrap, Equilibrium checks if it exists.
        let path = std::path::absolute(path).unwrap();

        let eq = Equilibrium::from_file(&path)?;

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

    /// Returns the `psip` coordinate data as a 1D array.
    pub fn psip_data(&self) -> Array1<f64> {
        Array1::from_vec(self.q_spline.xa.to_vec())
    }

    /// Returns the `q` data as a 1D array.
    pub fn q_data(&self) -> Array1<f64> {
        Array1::from_vec(self.q_spline.ya.to_vec())
    }

    /// Returns the `psi` data as a 1D array.
    pub fn psi_data(&self) -> Array1<f64> {
        Array1::from_vec(self.psi_spline.ya.to_vec())
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
        let path = PathBuf::from("./data.nc");
        Qfactor::from_dataset(&path, "akima").unwrap()
    }

    #[test]
    fn test_bfield_creation() {
        create_qfactor();
    }

    #[test]
    fn test_extraction_methods() {
        let q = create_qfactor();
        let _ = format!("{q:?}");

        assert_eq!(q.psip_data().shape(), [101]);
        assert_eq!(q.q_data().shape(), [101]);
        assert_eq!(q.psi_data().shape(), [101]);
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
