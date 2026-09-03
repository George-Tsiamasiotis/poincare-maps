//! Methods for data extraction from `netCDF` files and field names definitions.

/// The names of the netCDF fields (see [`convention`](https://dexter.tsiamasiotis.gr/netcdf)).
///
/// If the naming convention changes, this is the only file we must update.
pub mod netcdf_fields {

    // ================== Attributes ==================

    /// The netCDF's description.
    pub const NC_DESCRIPTION: &str = "description";
    /// The netCDF's creation date.
    pub const NC_DATE: &str = "date";
    /// The script used to create the netCDF file.
    pub const NC_SCRIPT: &str = "script";
    /// The netCDF's convention version ([`semver::Version`]).
    pub const NC_CONVENTION_VERSION: &str = "version";

    // ================== Scalars ==================

    /// The Magnetic field strength on the axis `B0` in **\[T\]**.
    pub const NC_BAXIS: &str = "baxis";
    /// The horizontal position of the magnetic axis `R0` in **\[m\]**.
    ///
    /// This field should be used for normalizations, rather than [`NC_RGEO`].
    pub const NC_RAXIS: &str = "raxis";
    /// The vertical position of the magnetic axis in **\[m\]**.
    pub const NC_ZAXIS: &str = "zaxis";
    /// The geometrical axis (device major radius) in **\[m\]**.
    pub const NC_RGEO: &str = "rgeo";

    // ================= Coordinates =================

    /// The boozer toroidal angle `θ` in **\[rads\]**.
    pub const NC_THETA: &str = "theta";
    /// The toroidal flux `ψ` in **Normalized Units**.
    pub const NC_PSI_NORM: &str = "psi_norm";
    /// The poloidal flux `ψp` in **Normalized Units**.
    pub const NC_PSIP_NORM: &str = "psip_norm";
    /// The radial coordinate `r` in **Normalized Units**.
    pub const NC_R_NORM: &str = "r_norm";
    /// The poloidal mode numbers `m`.
    pub const NC_M_MODES: &str = "m";
    /// The toroidal mode numbers `n`.
    pub const NC_N_MODES: &str = "n";
    /// The toroidal flux `ψ` in **\[Tm\]**.
    pub const NC_PSI: &str = "psi";
    /// The poloidal flux `ψp` in **\[Tm\]**.
    pub const NC_PSIP: &str = "psip";
    /// The radial coordinate `r` in **\[m\]**.
    pub const NC_R: &str = "r";

    // ================ 1D Variables ================

    /// The safety factor `q(ψ/ψp)`.
    pub const NC_Q: &str = "q";
    /// The covariant toroidal plasma current `g(ψ/ψp)` in **\[Tm\]**.
    pub const NC_G: &str = "g";
    /// The covariant poloidal plasma current `I(ψ/ψp)` in **\[Tm\]**.
    pub const NC_I: &str = "i";
    /// The covariant toroidal plasma current `g(ψ/ψp)` in **Normalized Units**.
    pub const NC_G_NORM: &str = "g_norm";
    /// The covariant poloidal plasma current `I(ψ/ψp)` in **Normalized Units**.
    pub const NC_I_NORM: &str = "i_norm";

    // ================ 2D Variables ================

    /// The magnetic field strength `B(ψ/ψp, θ)` in **\[T\]**.
    pub const NC_B: &str = "b";
    /// The magnetic field strength in `B(ψ/ψp, θ)` in **Normalized Units**.
    pub const NC_B_NORM: &str = "b_norm";
    /// The `R` coordinate with respect to boozer coordinates `R(ψ/ψp, θ)` in **\[m\]**.
    pub const NC_RLAB: &str = "rlab";
    /// The `Z` coordinate with respect to boozer coordinates `Z(ψ/ψp, θ)` in **\[m\]**.
    pub const NC_ZLAB: &str = "zlab";
    /// The VMEC to Boozer Jacobian `J(ψ/ψp, θ)` in **\[ m/T \]**.
    pub const NC_JACOBIAN: &str = "jacobian";

    // ================ 3D Variables ================

    /// The 3D array containing all the `α{m,n}(ψ/ψp)` 1D arrays in **Normalized Units**.
    pub const NC_ALPHAS_NORM: &str = "alphas_norm";
    /// The 3D array containing all the `α{m,n}(ψ/ψp)` 1D arrays in **\[m\]**.
    pub const NC_ALPHAS: &str = "alphas";
    /// The 3D array containing all the `φ{m,n}(ψ/ψp)` 1D arrays in **\[rads\]**.
    pub const NC_PHASES: &str = "phases";
}

// ===============================================================================================

use std::path::PathBuf;

use ndarray::{Array, Array1, Array2, Array3};
use netcdf::{Extents, NcTypeDescriptor, Variable};
use semver::Version;

use crate::NcError;

// ================== Testing netCDF files

/// Path to the netCDF file to be used for general testing.
/// Since this is used for testing the extraction methods, both fluxes should be good coordinates.
///
/// Must be *relative* to the **crate**'s `CARGO_MANIFEST_DIR`. Can be created with the
/// `tools/create_<>_test_netcdf.py`.
pub const TEST_NETCDF_PATH: &str = "test_netcdf.nc";

/// Test netcdf file with ψ the only good coordinate.
pub const TOROIDAL_TEST_NETCDF_PATH: &str = "toroidal_test_netcdf.nc";
/// Test netcdf file with ψp the only good coordinate.
pub const POLOIDAL_TEST_NETCDF_PATH: &str = "poloidal_test_netcdf.nc";

/// A `netCDF` file.
pub type NcFile = netcdf::File;

/// NetCDF-supported data types.
pub trait NcType: NcTypeDescriptor + Copy {}
impl NcType for f64 {}
impl NcType for f32 {}
impl NcType for i64 {}
impl NcType for i32 {}

/// Opens an [`NcFile`] from a given path.
///
/// # Example
/// ```
/// # use dexter_machine::*;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the path is not found or the `netCDF` file could not be opened.
pub fn open(path: &PathBuf) -> Result<NcFile, NcError> {
    if !path.exists() {
        return Err(NcError::FileNotFound(path.clone()));
    }
    match netcdf::open(path) {
        Ok(file) => Ok(file),
        Err(err) => Err(NcError::FileOpenError {
            path: path.clone(),
            err,
        }),
    }
}

/// Extracts an Attribute's value as a String.
///
/// # Example
/// ```
/// # use dexter_machine::extract::netcdf_fields;
/// # use dexter_machine::*;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
/// let date: String = extract::attribute(&file, netcdf_fields::NC_DATE)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Panics
///
/// Panics if [`NcFile::attribute`] panics.
///
/// # Errors
///
/// Returns an [`NcError`] if the [`netcdf::Attribute`] is not found or cannot be extracted.
pub fn attribute(file: &NcFile, name: &str) -> Result<String, NcError> {
    use netcdf::AttributeValue;

    let Some(attr) = file.attribute(name) else {
        return Err(NcError::AttributeNotFound(name.into()));
    };

    let Ok(value) = attr.value() else {
        return Err(NcError::AttributeValueError(name.into()));
    };

    // Even though all attributes are stored as strings, if there are any non-ascii characters,
    // they are stored as a vec containing a single string for some reason
    #[expect(clippy::wildcard_enum_match_arm, reason = "Attributes are strings")]
    match value {
        AttributeValue::Strs(items) => Ok(items.concat()),
        AttributeValue::Str(item) => Ok(item),
        _ => unreachable!(),
    }
}

/// Extracts the netCDF file's convention [`Semantic Version`].
///
/// # Example
/// ```
/// # use dexter_machine::*;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
/// let version: semver::Version = extract::version(&file)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`netcdf::Attribute`] is not found, or the version is not a valid
/// [`Semantic Version`].
///
/// [`Semantic Version`]: https://semver.org/
pub fn version(file: &NcFile) -> Result<Version, NcError> {
    let attr = attribute(file, netcdf_fields::NC_CONVENTION_VERSION)?;
    match Version::parse(&attr) {
        Ok(version) => Ok(version),
        Err(err) => Err(NcError::VersionError(err)),
    }
}

/// Checks if a [`Variable`] is not empty, e.g has a length of at least 1.
fn check_if_empty(var: &Variable) -> Result<(), NcError> {
    match var.len() {
        1.. => Ok(()),
        0 => Err(NcError::EmptyVariable(var.name())),
    }
}

/// Extracts a [`Variable`] named `name` from an [`NcFile`].
///
/// # Example
/// ```
/// # use dexter_machine::extract::netcdf_fields;
/// # use dexter_machine::*;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
///
/// let qfactor_var = extract::variable(&file, netcdf_fields::NC_Q)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found.
pub fn variable<'file>(file: &'file NcFile, name: &str) -> Result<Variable<'file>, NcError> {
    file.variable(name)
        .ok_or_else(|| NcError::VariableNotFound(name.into()))
}

// ===============================================================================================

/// Returns an [`Array<T, D>`] with the values of the [`Variable`] named `name`, in standard
/// layout.
fn array_nd<T, D>(file: &NcFile, name: &str) -> Result<Array<T, D>, NcError>
where
    T: NcType,
    D: ndarray::Dimension,
{
    let var = variable(file, name)?;
    check_if_empty(&var)?;

    let dyn_array = match var.get::<T, _>(Extents::All) {
        Ok(arr) => Ok(arr),
        Err(err) => Err(NcError::GetValues {
            name: var.name(),
            err,
        }),
    }?
    .into_dimensionality::<D>()?
    .as_standard_layout()
    .to_owned();

    Ok(dyn_array)
}

/// Extracts a scalar value of type `T` from an [`NcFile`].
///
/// # Example
/// ```
/// # use dexter_machine::extract::netcdf_fields;
/// # use dexter_machine::*;
/// # use ndarray::Array1;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
///
/// let baxis: f64 = extract::scalar(&file, netcdf_fields::NC_BAXIS)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn scalar<T: NcType>(file: &NcFile, name: &str) -> Result<T, NcError> {
    Ok(array_nd(file, name)?.into_scalar())
}

/// Extracts an [`Array1<T>`] value from an [`NcFile`].
///
/// The array is returned in [`standard layout`](method@ndarray::ArrayRef::as_standard_layout).
///
/// # Example
/// ```
/// # use dexter_machine::extract::netcdf_fields;
/// # use dexter_machine::*;
/// # use ndarray::Array1;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
///
/// let q_array: Array1<f64> = extract::array_1d(&file, netcdf_fields::NC_Q)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn array_1d<T: NcType>(file: &NcFile, name: &str) -> Result<Array1<T>, NcError> {
    array_nd(file, name)
}

/// Extracts an [`Array2<T>`] value from an [`NcFile`].
///
/// The array is returned in [`standard layout`](method@ndarray::ArrayRef::as_standard_layout).
///
/// # Example
/// ```
/// # use dexter_machine::extract::netcdf_fields;
/// # use dexter_machine::*;
/// # use ndarray::Array2;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
///
/// let b_norm_array: Array2<f64> = extract::array_2d(&file, netcdf_fields::NC_B_NORM)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn array_2d<T: NcType>(file: &NcFile, name: &str) -> Result<Array2<T>, NcError> {
    array_nd(file, name)
}

/// Extracts an [`Array3<T>`] value from an [`NcFile`].
///
/// The array is returned in [`standard layout`](method@ndarray::ArrayRef::as_standard_layout).
///
/// # Example
/// ```
/// # use dexter_machine::extract::netcdf_fields;
/// # use dexter_machine::*;
/// # use ndarray::Array3;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
///
/// let alphas_norm_array: Array3<f64> = extract::array_3d(&file, netcdf_fields::NC_ALPHAS_NORM)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the [`Variable`] is not found, is empty, or has a different shape.
pub fn array_3d<T: NcType>(file: &NcFile, name: &str) -> Result<Array3<T>, NcError> {
    array_nd(file, name)
}

/// Extracts the `α{m,n}(ψp)` and `φ{m,n}(ψp)` 1D arrays of the specified {m,n} mode.
///
/// # Example
/// ```
/// # use dexter_machine::*;
/// # use ndarray::Array1;
/// # use std::path::PathBuf;
/// #
/// let path = PathBuf::from("netcdf.nc");
/// let file = extract::open(&path)?;
///
/// let (mode32_alpha, mode32_phase) = extract::mode_arrays::<f64>(&file, 3, 2)?;
/// # Ok::<_, MachineError>(())
/// ```
///
/// # Errors
///
/// Returns an [`NcError`] if the `netCDF` file does not contain the {`m`, `n`} mode.
pub fn mode_arrays<T: NcType>(
    file: &NcFile,
    m_mode: i64,
    n_mode: i64,
) -> Result<(Array1<T>, Array1<T>), NcError> {
    let alpha_3d = array_3d::<T>(file, netcdf_fields::NC_ALPHAS_NORM)?;
    let phase_3d = array_3d::<T>(file, netcdf_fields::NC_PHASES)?;

    let m_index = get_logical_index(file, m_mode, netcdf_fields::NC_M_MODES)?;
    let n_index = get_logical_index(file, n_mode, netcdf_fields::NC_N_MODES)?;

    let alpha_1d = alpha_3d.slice(ndarray::s![m_index, n_index, ..]).to_owned();
    let phase_1d = phase_3d.slice(ndarray::s![m_index, n_index, ..]).to_owned();

    Ok((alpha_1d, phase_1d))
}

/// Returns the logical index of a mode's 1D arrays.
///
/// For example, if the `netCDF` file contains m = [-1, 0, 1, 2, 4], and we want the arrays
/// corresponding to m=1, we create the following index-mode mapping:
/// mmap = [
///     (0, -1),
///     (1, 0),
///     (2, 1),
///     (3, 2),
///     (4, 4),
/// ]
/// So for the mode m=1, we want the 2nd entry on the 3D array (`m_index` = 2).
///
/// It is agreed that the mode numbers are stored in an increasing order to avoid ambiguities.
#[expect(
    clippy::panic_in_result_fn,
    reason = "duplicate mode number in netCDF should be fatal"
)]
fn get_logical_index(file: &NcFile, mode: i64, field: &str) -> Result<usize, NcError> {
    let coord = array_1d::<i64>(file, field)?;
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
        .ok_or_else(|| NcError::FluteModeNotFound {
            which: field.to_lowercase(),
            mode,
        })
        .copied()
}

#[cfg(test)]
#[allow(unused_results)]
mod test {
    use super::netcdf_fields::*;
    use super::*;

    fn open_test_file() -> NcFile {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        open(&path).unwrap()
    }

    #[test]
    fn netcdf_open() {
        open_test_file();
    }

    #[test]
    fn netcdf_all_attributes_extraction() {
        let file = open_test_file();

        attribute(&file, NC_DESCRIPTION).unwrap();
        attribute(&file, NC_DATE).unwrap();
        attribute(&file, NC_SCRIPT).unwrap();
        attribute(&file, NC_CONVENTION_VERSION).unwrap();
        version(&file).unwrap();
    }

    #[test]
    fn extract_raw_variable() {
        let file = open_test_file();

        let var = variable(&file, NC_BAXIS).unwrap();
        let var_units = var.attribute("units").unwrap();
        assert_eq!(var_units.name(), "units");

        use netcdf::AttributeValue;
        assert_eq!(
            var_units.value().unwrap(),
            AttributeValue::Str("[T]".into())
        );
    }

    #[test]
    fn netcdf_all_scalars_extraction() {
        let file = open_test_file();

        scalar::<f64>(&file, NC_BAXIS).unwrap();
        scalar::<f64>(&file, NC_RAXIS).unwrap();
        scalar::<f64>(&file, NC_ZAXIS).unwrap();
        scalar::<f64>(&file, NC_RGEO).unwrap();
    }

    #[test]
    fn netcdf_all_1d_arrays_extraction() {
        let file = open_test_file();

        array_1d::<f64>(&file, NC_THETA).unwrap();
        array_1d::<f64>(&file, NC_PSIP_NORM).unwrap();
        array_1d::<f64>(&file, NC_PSI_NORM).unwrap();
        array_1d::<f64>(&file, NC_R_NORM).unwrap();
        array_1d::<i64>(&file, NC_M_MODES).unwrap();
        array_1d::<i64>(&file, NC_N_MODES).unwrap();
        array_1d::<f64>(&file, NC_PSIP).unwrap();
        array_1d::<f64>(&file, NC_PSI).unwrap();
        array_1d::<f64>(&file, NC_R).unwrap();

        array_1d::<f64>(&file, NC_Q).unwrap();
        array_1d::<f64>(&file, NC_G).unwrap();
        array_1d::<f64>(&file, NC_I).unwrap();
        array_1d::<f64>(&file, NC_G_NORM).unwrap();
        array_1d::<f64>(&file, NC_I_NORM).unwrap();
    }

    #[test]
    fn netcdf_all_2d_arrays_extraction() {
        let file = open_test_file();

        array_2d::<f64>(&file, NC_B).unwrap();
        array_2d::<f64>(&file, NC_B_NORM).unwrap();
        array_2d::<f64>(&file, NC_RLAB).unwrap();
        array_2d::<f64>(&file, NC_ZLAB).unwrap();
        array_2d::<f64>(&file, NC_JACOBIAN).unwrap();
    }

    #[test]
    fn netcdf_all_3d_arrays_extraction() {
        let file = open_test_file();

        array_3d::<f64>(&file, NC_ALPHAS_NORM).unwrap();
        array_3d::<f64>(&file, NC_ALPHAS).unwrap();
        array_3d::<f64>(&file, NC_PHASES).unwrap();
    }

    /// WARN: Make sure this test is up to date with the stub netcdf file.
    /// We inspect the (2,2) mode, which corresponds to the indices (0, 1).
    #[test]
    fn nc_flute_mode_extraction_values() {
        let file = open_test_file();

        let alpha_3d = array_3d::<f64>(&file, NC_ALPHAS_NORM).unwrap();
        let phase_3d = array_3d::<f64>(&file, NC_PHASES).unwrap();

        // Cast to i64 to avoid float comparisons
        use ndarray::s;
        // Index by logical index
        assert_eq!(
            alpha_3d.slice(s![0, 1, 0]).into_scalar().to_owned() as i64,
            1111,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            phase_3d.slice(s![0, 1, 0]).into_scalar().to_owned() as i64,
            9999,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            alpha_3d.slice(s![0, 1, -1]).into_scalar().to_owned() as i64,
            11111,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            phase_3d.slice(s![0, 1, -1]).into_scalar().to_owned() as i64,
            99999,
            "Is this test up to date with the stub netcdf file?"
        );

        // Index by mode number
        assert_eq!(
            mode_arrays::<f64>(&file, 2, 2)
                .unwrap()
                .0
                .first()
                .unwrap()
                .to_owned() as i64,
            1111,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            mode_arrays::<f64>(&file, 2, 2)
                .unwrap()
                .0
                .last()
                .unwrap()
                .to_owned() as i64,
            11111,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            mode_arrays::<f64>(&file, 2, 2)
                .unwrap()
                .1
                .first()
                .unwrap()
                .to_owned() as i64,
            9999,
            "Is this test up to date with the stub netcdf file?"
        );
        assert_eq!(
            mode_arrays::<f64>(&file, 2, 2)
                .unwrap()
                .1
                .last()
                .unwrap()
                .to_owned() as i64,
            99999,
            "Is this test up to date with the stub netcdf file?"
        );
    }

    #[test]
    fn netcdf_errors() {
        let file = open_test_file();

        assert!(matches!(
            dbg!(open(&PathBuf::from("not a path"))),
            Err(NcError::FileNotFound(..))
        ));

        assert!(matches!(
            dbg!(open(&PathBuf::from("/tmp"))),
            Err(NcError::FileOpenError { .. })
        ));

        assert!(matches!(
            dbg!(attribute(&file, "not an attribute")),
            Err(NcError::AttributeNotFound(..))
        ));

        assert!(matches!(
            dbg!(variable(&file, "not an attribute")),
            Err(NcError::VariableNotFound(..))
        ));

        assert!(matches!(
            dbg!(scalar::<f64>(&file, "not a name")),
            Err(NcError::VariableNotFound(..))
        ));

        assert!(matches!(
            dbg!(dbg!(array_1d::<f64>(&file, NC_B_NORM))),
            Err(NcError::ShapeError { .. })
        ));

        assert!(matches!(
            dbg!(mode_arrays::<f64>(&file, 1000, -2000)),
            Err(NcError::FluteModeNotFound { .. })
        ));
    }
}
