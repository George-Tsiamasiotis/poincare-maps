use std::f64::consts::TAU;
use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline, make_spline};

use crate::Result;

use ndarray::Array1;
use numpy::{PyArray1, ToPyArray};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// Single perturbation harmonic reconstructed from a netCDF file.
///
/// The harmonic has the form of `α(ψp) * cos(mθ-nζ+φ0)`, where `α(ψp)` is calculated by
/// interpolation over some numerical data.
#[pyclass(frozen, immutable_type)]
pub struct Harmonic {
    /// Path to the netCDF file.
    #[pyo3(get)]
    pub path: PathBuf,
    /// Interpolation type.
    #[pyo3(get)]
    pub typ: String,

    /// Spline over the perturbation amplitude `α` data, as a function of ψp.
    pub a_spline: DynSpline<f64>,
    /// The `θ` frequency number.
    pub m: f64,
    /// The `ζ` frequency number.
    pub n: f64,
    /// The initial phase of the harmonic.
    pub phase: f64,

    /// The maximum value of the `α` values.
    #[pyo3(get)]
    pub amax: f64,
    /// The value of the poloidal angle ψp at the wall.
    #[pyo3(get)]
    pub psip_wall: f64,
}

/// Wrapper methods exposed to Python.
#[pymethods]
impl Harmonic {
    /// Creates a new [`Harmonic`]
    ///
    /// Wrapper around [`Harmonic::from_dataset`]. This is a workaround to return a [`PyErr`].
    #[coverage(off)]
    #[new]
    pub fn new_py(path: &str, typ: &str, m: f64, n: f64, phase: f64) -> PyResult<Self> {
        let path = PathBuf::from(path);
        match Self::from_dataset(&path, typ, m, n, phase) {
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

    /// Returns the `α` amplitude data as a Numpy 1D array.
    #[coverage(off)]
    #[pyo3(name = "a_data")]
    pub fn a_data_py<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.a_data().to_pyarray(py)
    }

    #[coverage(off)]
    pub fn __repr__(&self) -> String {
        format!("{:#?}", &self)
    }
}

impl Harmonic {
    /// Constructs a [`Harmonic`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// The spline is only over the amplitude `α`, of the perturbation, and the rest of the
    /// exrpession is analytic.
    ///
    /// # Note
    ///
    /// The value `ψ = 0.0` is prepended at the ψ data array, and the first value of the q array is
    /// prepended (duplicated) in the q array, to assure correct interpolation near the magnetic axis.
    ///
    /// # Example
    /// ```
    /// # use poincare_maps::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("./data.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str, m: f64, n: f64, phase: f64) -> Result<Self> {
        use rsl_interpolation::*;
        use tokamak_netcdf::variable_names::*;
        use tokamak_netcdf::*;

        // Make path absolute. Just unwrap, Equilibrium checks if it exists.
        let path = std::path::absolute(path).unwrap();

        let eq = Equilibrium::from_file(&path)?;

        // Add 0.0 manualy, which corresponds to α(0).
        let psip_data = extract_var_with_axis_value(&eq.file, PSIP_COORD, 0.0)?
            .as_standard_layout()
            .to_vec();
        let psip_wall = psip_data.last().copied().unwrap();

        // TODO: update
        let mut a_data = Array1::zeros(psip_data.len());
        for i in 0..psip_data.len() {
            a_data[i] = gaussian(psip_data[i], psip_wall)
        }
        let amax = *a_data
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&0.0);

        let a_spline = make_spline(typ, &psip_data, &a_data.to_vec())?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            a_spline,
            psip_wall,
            m,
            n,
            phase: phase % TAU,
            amax,
        })
    }

    /// Returns the `psip` coordinate data as a 1D array.
    pub fn psip_data(&self) -> Array1<f64> {
        Array1::from_vec(self.a_spline.xa.to_vec())
    }

    /// Returns the `α` perturbation amplitude data as a 1D array.
    pub fn a_data(&self) -> Array1<f64> {
        Array1::from_vec(self.a_spline.ya.to_vec())
    }
}

impl Harmonic {
    /// Calculates the harmonic `a(ψp) * <analytical term>`.
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
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let h = harmonic.h(0.015, 2.0*PI, 0.0, &mut psi_acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn h(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let term = (self.m * theta - self.n * zeta + self.phase).cos();
        Ok(self.a_spline.eval(psip, acc)? * term)
    }

    /// Calculates the harmonic derivative `𝜕h/𝜕ψp`.
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
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dh_dpsip = harmonic.dh_dpsip(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dh_dpsip(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let term = (self.m * theta - self.n * zeta + self.phase).cos();
        Ok(self.a_spline.eval_deriv(psip, acc)? * term)
    }

    /// Calculates the harmonic derivative `𝜕h/𝜕θ`.
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
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dh_dtheta = harmonic.dh_dtheta(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dh_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        let term = -self.m * (self.m * theta - self.n * zeta + self.phase).sin();
        Ok(self.a_spline.eval(psip, acc)? * term)
    }

    /// Calculates the perturbation derivative `𝜕h/𝜕ζ`.
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
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dh_dzeta = harmonic.dh_dzeta(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dh_dzeta(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        let term = self.n * (self.m * theta - self.n * zeta + self.phase).sin();
        Ok(self.a_spline.eval(psip, acc)? * term)
    }

    /// Calculates the perturbation derivative `𝜕h/𝜕t`.
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
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dh_dt = harmonic.dh_dt(0.015, 2.0*PI, PI, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    #[allow(unused_variables)]
    pub fn dh_dt(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }
}

/// A simple gaussian distribution to emulate reconstructed perturbations.
/// TODO: remove
fn gaussian(psip: f64, psip_wall: f64) -> f64 {
    use std::f64::consts::TAU;

    let scale = 1e-4;
    let mu = psip_wall / 2.0;
    let sigma = psip_wall / 4.0;

    scale * (1.0 / (TAU * sigma).sqrt()) * (-(psip - mu).powi(2) / (2.0 * sigma.powi(2))).exp()
}

impl Clone for Harmonic {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            typ: self.typ.clone(),
            a_spline: make_spline(&self.typ, &self.a_spline.xa, &self.a_spline.ya)
                .expect("Could not clone spline."),
            m: self.m,
            n: self.n,
            phase: self.phase,
            amax: self.amax,
            psip_wall: self.psip_wall,
        }
    }
}

impl std::fmt::Debug for Harmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Harmonic")
            .field("path", &self.path)
            .field("typ", &self.typ)
            .field("ψp_wall", &format!("{:.7}", self.psip_wall))
            .field("m", &self.m)
            .field("n", &self.n)
            .field("φ", &self.phase)
            .field("α_max", &self.amax)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_data_extraction() {
        let path = PathBuf::from("./data.nc");
        let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0).unwrap();

        assert_eq!(harmonic.psip_data().ndim(), 1);
        assert_eq!(harmonic.a_data().ndim(), 1);
    }

    #[test]
    fn test_harmonic_misc() {
        let path = PathBuf::from("./data.nc");
        let harmonic = Harmonic::from_dataset(&path, "akima", 3.0, 2.0, 0.0).unwrap();
        let _ = harmonic.clone();
        let _ = format!("{harmonic:?}");
    }
}
