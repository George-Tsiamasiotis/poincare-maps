//! Common getter method implementations for machine objects.

/// Creates a getter method for extracting the flat Vec data as an Array2.
/// The Vec is assumed to be in Fortran order, since it is intended for use by the splines.
#[doc(hidden)]
#[macro_export]
macro_rules! fortran_vec_to_carray2d_impl {
    ($meth_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values as a 2D array." ]
        #[must_use]
        pub fn $meth_name(&self) -> Array2<f64> {
            // Array is in Fortran order, so we must reverse the shape
            let actual_shape = self.shape();
            let shape = (actual_shape.1, actual_shape.0);
            Array2::from_shape_vec(shape, self.$($field).+.clone())
                .expect("Shape is correct by definition")
                .reversed_axes()
        }
    }
}

/// Generates getters for the fluxes' values.
#[doc(hidden)]
#[macro_export]
macro_rules! fluxes_values_array_getter_impl {
    () => {
        /// Returns the toroidal flux's values as a 1D array, if they exist.
        #[must_use]
        pub fn psi_array(&self) -> Option<Array1<f64>> {
            self.psi
                .values()
                .map(|values| Array1::from(Vec::from(values)))
        }

        /// Returns the poloidal flux's values as a 1D array, if they exist.
        #[must_use]
        pub fn psip_array(&self) -> Option<Array1<f64>> {
            self.psip
                .values()
                .map(|values| Array1::from(Vec::from(values)))
        }
    };
}

/// Generates getters for the last closed flux surfaces, for `Nc` types that do not
/// implement `Qfactor`.
#[doc(hidden)]
#[macro_export]
macro_rules! lcfs_getter_impl {
    () => {
        /// Returns the value of the last closed toroidal flux `ψ_last`.
        #[must_use]
        pub fn psi_last(&self) -> Option<f64> {
            self.psi.last_value()
        }

        /// Returns the value of the last closed poloidal flux `ψp_last`.
        #[must_use]
        pub fn psip_last(&self) -> Option<f64> {
            self.psip.last_value()
        }
    };
}

/// Generates getters for a [`crate::ModeCache`] implementor's hits and misses
#[doc(hidden)]
#[macro_export]
macro_rules! mode_cache_getters_impl {
    ($obj: ident) => {
        fn hits(&self) -> usize {
            self.hits
        }

        fn misses(&self) -> usize {
            self.misses
        }

        fn cache(&mut self) -> &mut [f64] {
            &mut self.cache
        }

        fn params(&mut self) -> &[f64] {
            &self.params
        }
    };
}

/// Generates a getter for the object's path to netCDF file.
#[doc(hidden)]
#[macro_export]
macro_rules! netcdf_path_getter_impl {
    () => {
        /// Returns the netCDF file's path.
        #[must_use]
        pub fn path(&self) -> PathBuf {
            self.path.clone()
        }
    };
}

/// Generates a getter for the object's Interpolation types
#[doc(hidden)]
#[macro_export]
macro_rules! interp_type_getter_impl {
    // 1D Interpolation
    (1) => {
        /// Returns the interpolation type.
        #[must_use]
        pub fn interp_type(&self) -> Interpolation1dType {
            self.interp_type
        }
    };
    // 2D Interpolation
    (2) => {
        /// Returns the interpolation type.
        #[must_use]
        pub fn interp_type(&self) -> Interpolation2dType {
            self.interp_type
        }
    };
}

/// Generates a getter for the object's netCDF convention version.
#[doc(hidden)]
#[macro_export]
macro_rules! netcdf_version_getter_impl {
    () => {
        /// Returns the object's [`Version`](semver::Version).
        #[must_use]
        pub fn netcdf_version(&self) -> semver::Version {
            self.netcdf_version.clone()
        }
    };
}
