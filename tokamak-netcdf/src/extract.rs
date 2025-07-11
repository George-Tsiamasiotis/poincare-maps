//! Functions for extracting and checking data from the NetCDF file.

use crate::{NcError, Result};
use ndarray::{Array1, Array2, ArrayView, Axis, array};

/// Extracts a `Variable` fron a NetCDF file.
fn extract_variable<'a>(f: &'a netcdf::File, name: &'a str) -> Result<netcdf::Variable<'a>> {
    f.variable(name)
        .ok_or(NcError::VariableNotFound(name.into()))
}

/// Checks if a `Variable` is empty.
fn check_if_empty(var: &netcdf::Variable) -> Result<()> {
    match var.len() {
        1.. => Ok(()),
        0 => Err(NcError::EmptyVariable(var.name().into())),
    }
}

/// Extracts a scalar (0D) `Variable`'s value.
pub(crate) fn extract_scalar<T>(f: &netcdf::File, name: &str) -> Result<T>
where
    T: netcdf::NcTypeDescriptor + std::marker::Copy,
{
    use crate::NcError::*;

    let var = extract_variable(f, name)?;
    check_if_empty(&var)?;

    // `var.dimensions()` is () for netcdf's scalar `Variables`. This is probably equivalent to
    // `var.len() == 0`
    if !var.dimensions().is_empty() {
        return Err(NotScalar(name.into()));
    }

    match var.get_value::<T, _>(..) {
        Ok(value) => Ok(value),
        Err(err) => Err(NcError::GetValuesError {
            name: var.name().into(),
            source: err,
        }),
    }
}

/// Extracts a 1D `Variable` and returns its values.
pub(crate) fn extract_1d_var<T>(f: &netcdf::File, name: &str) -> Result<Array1<T>>
where
    T: netcdf::NcTypeDescriptor + std::marker::Copy + std::default::Default,
{
    let var = extract_variable(f, name)?;
    check_if_empty(&var)?;

    if var.dimensions().len() != 1 {
        return Err(NcError::Not1D(var.name().into()));
    }

    let mut data = Array1::<T>::default(var.len());

    match var.get_into(data.view_mut(), ..) {
        Ok(()) => Ok(data),
        Err(err) => Err(NcError::GetValuesError {
            name: var.name().into(),
            source: err,
        }),
    }
}

/// Extracts a 2D `Variable` and returns its values as an `ndarray`.
pub(crate) fn extract_2d_var<T>(f: &netcdf::File, name: &str) -> Result<Array2<T>>
where
    T: netcdf::NcTypeDescriptor + std::marker::Copy + std::default::Default,
{
    let var = extract_variable(f, name)?;
    check_if_empty(&var)?;

    if var.dimensions().len() != 2 {
        return Err(NcError::Not2D(var.name().into()));
    }

    // Dimension order is (ψ,θ).
    let dims = var.dimensions().to_vec();
    let shape = (dims[0].len(), dims[1].len());

    let mut data = Array2::<T>::default(shape);

    match var.get_into(data.view_mut(), (.., ..)) {
        Ok(()) => Ok(data),
        Err(err) => Err(NcError::GetValuesError {
            name: var.name().into(),
            source: err,
        }),
    }
}

/// Extracts a variable from the NetCDF file and prepends the first value (value closest to the
/// magnetic axis) at index 0.
pub(crate) fn extract_var_with_first_axis_value<T>(
    f: &netcdf::File,
    name: &str,
) -> Result<Array1<T>>
where
    T: netcdf::NcTypeDescriptor + std::marker::Copy + std::default::Default,
{
    let arr: Array1<T> = extract_1d_var(f, name)?;
    extract_var_with_axis_value(f, name, arr[0])
}

/// Extracts a variable from the NetCDF file and prepends `element` at index 0.
pub(crate) fn extract_var_with_axis_value<T>(
    f: &netcdf::File,
    name: &str,
    element: T,
) -> Result<Array1<T>>
where
    T: netcdf::NcTypeDescriptor + std::marker::Copy + std::default::Default,
{
    let arr: Array1<T> = extract_1d_var(f, name)?;
    let view = ArrayView::from(&arr);
    let mut prepend: Array1<T> = array![element];
    // This is not expected to fail since both arrays are guranteed to be of the same shape (1,).
    match prepend.append(Axis(0), view) {
        Ok(()) => Ok(prepend),
        Err(_) => unreachable!("Shape mismatch in prepending axis value."),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use NcError::*;

    static VAR_LENGTH: usize = 5;

    /// Creates a phony NetCDF file for use across the tests.
    fn phony_netcdf() -> std::result::Result<netcdf::FileMut, netcdf::Error> {
        let path = std::env::temp_dir().join("phony.nc");
        let path_str = path.to_str().unwrap();

        let mut f = netcdf::create(path_str)?;
        std::fs::remove_file(path).unwrap();

        f.add_dimension("dim1", VAR_LENGTH)?;
        f.add_dimension("dim2", VAR_LENGTH)?;
        f.add_variable::<f64>("var", &["dim2"])?;

        f.add_dimension("empty_dim", 0)?;
        f.add_variable::<f64>("empty_var", &["empty_dim"])?;

        f.add_variable::<f64>("2dvar", &["dim1", "dim2"])?;
        f.add_variable::<i32>("int_var", &["dim1"])?;

        f.add_variable::<i32>("number", &[])?
            .put_values(&[18], ..)?;
        Ok(f)
    }

    #[test]
    fn test_extract_variable() {
        let f = phony_netcdf().unwrap();
        extract_variable(&f, "var").unwrap();
        assert!(matches!(
            extract_variable(&f, "not_a_var").unwrap_err(),
            VariableNotFound(_)
        ));
    }

    #[test]
    fn test_check_if_empty() -> Result<()> {
        let f = phony_netcdf().unwrap();
        let var = extract_variable(&f, "var")?;
        let empty_var = extract_variable(&f, "empty_var")?;

        assert_eq!(var.len(), VAR_LENGTH);
        assert_eq!(empty_var.len(), 0);
        assert!(matches!(
            check_if_empty(&empty_var).unwrap_err(),
            EmptyVariable(_)
        ));
        Ok(())
    }

    #[test]
    fn test_extract_scalar() -> Result<()> {
        let f = phony_netcdf().unwrap();
        let scalar: i32 = extract_scalar(&f, "number")?;
        let not_a_scalar = extract_scalar::<f64>(&f, "var");

        assert_eq!(scalar, 18i32);
        assert!(matches!(not_a_scalar.unwrap_err(), NotScalar(_)));
        Ok(())
    }

    #[test]
    fn test_extract_1d_var() {
        let f = phony_netcdf().unwrap();
        let values1d = extract_1d_var::<f64>(&f, "var");
        let values2d = extract_1d_var::<f64>(&f, "2dvar");
        let empty_values = extract_1d_var::<f64>(&f, "empty_var");
        let err_values = extract_1d_var::<f64>(&f, "not_a_var");

        values1d.unwrap();
        assert!(matches!(values2d.unwrap_err(), Not1D(_)));
        assert!(matches!(empty_values.unwrap_err(), EmptyVariable(_)));
        assert!(matches!(err_values.unwrap_err(), VariableNotFound(_)));
    }

    #[test]
    fn test_ectract_2d_var() {
        let f = phony_netcdf().unwrap();
        let values2d = extract_2d_var::<f64>(&f, "2dvar");
        let values1d = extract_2d_var::<f64>(&f, "var");
        let empty_values = extract_2d_var::<f64>(&f, "empty_var");
        let err_values = extract_2d_var::<f64>(&f, "not_a_var");

        values2d.unwrap();
        assert!(matches!(values1d.unwrap_err(), Not2D(_)));
        assert!(matches!(empty_values.unwrap_err(), EmptyVariable(_)));
        assert!(matches!(err_values.unwrap_err(), VariableNotFound(_)));
    }

    #[test]
    fn test_axis_value() -> Result<()> {
        let mut f = phony_netcdf().unwrap();
        let data: [i32; VAR_LENGTH] = [2, 3, 4, 5, 6];

        f.variable_mut("int_var")
            .expect("Error extracting mutable variable.")
            .put_values(&data, ..)
            .expect("Error putting values to variable");

        assert_eq!(
            Array1::<i32>::from_vec(vec![1, 2, 3, 4, 5, 6]),
            extract_var_with_axis_value(&f, "int_var", 1)?
        );
        assert_eq!(
            Array1::<i32>::from_vec(vec![2, 2, 3, 4, 5, 6]),
            extract_var_with_first_axis_value(&f, "int_var")?
        );
        Ok(())
    }
}
