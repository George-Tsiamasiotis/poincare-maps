//! Methods for data extraction from NetCDF file.
//!
//! - For spline creation, only the `extract_<>d_array()` functions are necessary.
//!
//! - `extract_variable` provides access to the underlying [`Variable`] and it's [`netcdf::Attribute`]s.

use std::path::PathBuf;

use config::netcdf_fields::{ALPHAS, M, N, PHASES};
use ndarray::{Array, Array1, Array2, Array3};
use netcdf::{Extents, File, NcTypeDescriptor, Variable};

use crate::NcError;

type Result<T> = std::result::Result<T, NcError>;

/// NetCDF-supported data types.
pub trait NcType: NcTypeDescriptor + Copy {}
impl NcType for f64 {}
impl NcType for f32 {}
impl NcType for i64 {}
impl NcType for i32 {}

/// Opens a [`File`] from a given path.
///
/// # Example
/// ```
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = open(&path)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the path is not found or the NetCDF could not be opened.
pub fn open(path: &PathBuf) -> Result<File> {
    if !path.exists() {
        return Err(NcError::FileNotFound(path.clone()));
    }
    match netcdf::open(path) {
        Ok(f) => Ok(f),
        Err(err) => Err(NcError::FileOpenError {
            path: path.clone(),
            err,
        }),
    }
}

/// Checks if a [`Variable`] is not empty, e.g has a length of at least 1.
fn check_if_empty(var: &Variable) -> Result<()> {
    match var.len() {
        1.. => Ok(()),
        0 => Err(NcError::EmptyVariable(var.name())),
    }
}

/// Extracts a [`Variable`] named `name` from a [`File`].
///
/// # Example
/// ```
/// # use config::netcdf_fields::*;
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use netcdf::Variable;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = extract::open(&path)?;
///
/// let qfactor_var: Variable = extract_variable(&f, BAXIS)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found.
pub fn extract_variable<'f>(f: &'f File, name: &str) -> Result<Variable<'f>> {
    f.variable(name)
        .ok_or(NcError::VariableNotFound(name.into()))
}

// ===============================================================================================

/// Returns an [`Array<T, D>`] with the values of the [`Variable`] named `name`.
fn extract_array<T, D>(f: &File, name: &str) -> Result<Array<T, D>>
where
    T: NcType,
    D: ndarray::Dimension,
{
    let var = extract_variable(f, name)?;
    check_if_empty(&var)?;

    let dyn_array = match var.get::<T, _>(Extents::All) {
        Ok(arr) => Ok(arr),
        Err(err) => Err(NcError::GetValues {
            name: var.name(),
            err,
        }),
    }?
    .into_dimensionality::<D>()?;

    Ok(dyn_array)
}

/// Extracts a scalar value of type `T` from a [`File`].
///
/// # Example
/// ```
/// # use config::netcdf_fields::*;
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use ndarray::Array1;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = extract::open(&path)?;
///
/// let baxis: f64 = extract_scalar(&f, BAXIS)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn extract_scalar<T: NcType>(f: &File, name: &str) -> Result<T> {
    Ok(extract_array(f, name)?.into_scalar())
}

/// Extracts an [`Array1<T>`] value from a [`File`].
///
/// # Example
/// ```
/// # use config::netcdf_fields::*;
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use ndarray::Array1;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = extract::open(&path)?;
///
/// let qfactor: Array1<f64> = extract_1d_array(&f, Q)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn extract_1d_array<T: NcType>(f: &File, name: &str) -> Result<Array1<T>> {
    extract_array(f, name)
}

/// Extracts an [`Array2<T>`] value from a [`File`].
///
/// # Example
/// ```
/// # use config::netcdf_fields::*;
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use ndarray::Array2;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = extract::open(&path)?;
///
/// let b_norm: Array2<f64> = extract_2d_array(&f, B_NORM)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn extract_2d_array<T: NcType>(f: &File, name: &str) -> Result<Array2<T>> {
    extract_array(f, name)
}

/// Extracts an [`Array3<T>`] value from a [`File`].
///
/// # Example
/// ```
/// # use config::netcdf_fields::*;
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use ndarray::Array3;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = extract::open(&path)?;
///
/// let alphas_norm: Array3<f64> = extract_3d_array(&f, ALPHAS_NORM)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn extract_3d_array<T: NcType>(f: &File, name: &str) -> Result<Array3<T>> {
    extract_array(f, name)
}

/// Extracts the `α{m,n}(ψp)` and `φ{m,n}(ψp)` 1D arrays of the specified {m,n} mode.
///
/// # Example
/// ```
/// # use config::netcdf_fields::*;
/// # use equilibrium::extract::*;
/// # use equilibrium::*;
/// # use ndarray::Array1;
/// # use std::path::PathBuf;
/// #
/// # fn main() -> Result<()> {
/// let path = PathBuf::from("../data/stub_netcdf.nc");
/// let f = extract::open(&path)?;
///
/// let (harmonic32_alpha, harmonic32_phase): (Array1<f64>, Array1<f64>) =
///     extract_harmonic_arrays(&f, 3, 2)?;
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the NetCDF file does not contain the {`m`, `n`} harmonic.
pub fn extract_harmonic_arrays<T: NcType>(
    f: &File,
    m: i64,
    n: i64,
) -> Result<(Array1<T>, Array1<T>)> {
    let alpha_3d = extract_3d_array::<T>(f, ALPHAS)?;
    let phase_3d = extract_3d_array::<T>(f, PHASES)?;

    let m_index = get_logical_index(f, m, M)?;
    let n_index = get_logical_index(f, n, N)?;

    let alpha_1d = alpha_3d.slice(ndarray::s![m_index, n_index, ..]).to_owned();
    let phase_1d = phase_3d.slice(ndarray::s![m_index, n_index, ..]).to_owned();

    Ok((alpha_1d, phase_1d))
}

/// Returns the logical index of a harmonic's 1D arrays.
///
/// For example, if the NetCDF file contains m = [-1, 0, 1, 2, 4], and we want the arrays
/// corresponding to m=1, we create the following index-mode mapping:
/// mmap = [
///     (0, -1),
///     (1, 0),
///     (2, 1),
///     (3, 2),
///     (4, 4),
/// ]
/// So for the mode m=1, we want the 2nd entry on the 3D array (m_index = 2).
fn get_logical_index(f: &File, mode: i64, field: &str) -> Result<usize> {
    let coord = extract_1d_array::<i64>(f, field)?;
    let map: Vec<(usize, &i64)> = coord.indexed_iter().collect(); // create mapping
    let pair = map
        .iter()
        .filter(|tuple| *tuple.1 == mode) // pick the (index, mode) entry we want
        .map(|tuple| tuple.0) // drop index
        .collect::<Vec<usize>>();
    // There should be at most 1 entry at this point, otherwise there is something
    // seriously wrong with the data.
    assert!(pair.len() <= 1, "Duplicate mode numbers found");

    pair.first()
        .ok_or(NcError::HarmonicModeNotFound {
            which: "field".to_lowercase(),
            mode,
        })
        .copied()
}

#[cfg(test)]
mod test {
    use super::*;
    use config::STUB_NETCDF_PATH;
    use config::netcdf_fields::*;

    fn open_test_file() -> netcdf::File {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        open(&path).unwrap()
    }

    #[test]
    fn test_netcdf_open() {
        open_test_file();
    }

    #[test]
    fn test_netcdf_all_scalars_extraction() {
        let f = open_test_file();

        extract_scalar::<f64>(&f, BAXIS).unwrap();
        extract_scalar::<f64>(&f, RAXIS).unwrap();
    }

    #[test]
    fn test_netcdf_all_1d_arrays_extraction() {
        let f = open_test_file();

        extract_1d_array::<f64>(&f, THETA).unwrap();
        extract_1d_array::<f64>(&f, PSIP).unwrap();
        extract_1d_array::<f64>(&f, PSI).unwrap();
        extract_1d_array::<i64>(&f, M).unwrap();
        extract_1d_array::<i64>(&f, N).unwrap();

        extract_1d_array::<f64>(&f, Q).unwrap();
        extract_1d_array::<f64>(&f, G).unwrap();
        extract_1d_array::<f64>(&f, I).unwrap();
        extract_1d_array::<f64>(&f, G_NORM).unwrap();
        extract_1d_array::<f64>(&f, I_NORM).unwrap();
    }

    #[test]
    fn test_netcdf_all_2d_arrays_extraction() {
        let f = open_test_file();

        extract_2d_array::<f64>(&f, B).unwrap();
        extract_2d_array::<f64>(&f, B_NORM).unwrap();
        extract_2d_array::<f64>(&f, R).unwrap();
        extract_2d_array::<f64>(&f, Z).unwrap();
    }

    #[test]
    fn test_netcdf_all_3d_arrays_extraction() {
        let f = open_test_file();

        extract_3d_array::<f64>(&f, ALPHAS_NORM).unwrap();
        extract_3d_array::<f64>(&f, ALPHAS).unwrap();
        extract_3d_array::<f64>(&f, PHASES).unwrap();
    }

    /// WARN: Make sure this test is up to date with the stub netcdf file.
    #[test]
    fn test_netcdf_harmonic_extraction_values() {
        let f = open_test_file();

        let alpha_3d = extract_3d_array::<f64>(&f, ALPHAS).unwrap();
        let phase_3d = extract_3d_array::<f64>(&f, PHASES).unwrap();

        // Cast to i64 to avoid float comparisons
        use ndarray::s;
        assert_eq!(
            alpha_3d.slice(s![2, 3, 0]).into_scalar().to_owned() as i64,
            1111,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            phase_3d.slice(s![2, 3, 0]).into_scalar().to_owned() as i64,
            9999,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            alpha_3d.slice(s![2, 3, -1]).into_scalar().to_owned() as i64,
            11111,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            phase_3d.slice(s![2, 3, -1]).into_scalar().to_owned() as i64,
            99999,
            "Is this test up to date with the stub netcdf file?"
        );
    }

    #[test]
    fn test_netcdf_errors() {
        let f = open_test_file();

        assert!(matches!(
            open(&PathBuf::from("not a path")),
            Err(NcError::FileNotFound(..))
        ));

        assert!(matches!(
            extract_scalar::<f64>(&f, "not a name"),
            Err(NcError::VariableNotFound(..))
        ));

        assert!(matches!(
            extract_harmonic_arrays::<f64>(&f, 1000, -2000),
            Err(NcError::HarmonicModeNotFound { .. })
        ));
    }
}
