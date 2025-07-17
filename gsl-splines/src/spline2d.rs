use crate::{Accelerator, Interpolation2dType, RgslSpline2d, SplineError};
use crate::{Result, RgslInterp2dType};
use ndarray::{Array1, Array2};
use std::cell::RefCell;
use std::fmt::Display;
use std::rc::Rc;

/// A 2D Spline of a specific [`type`] that stores the data points and can share both [`Accelerators`] with
/// other 1D or 2D splines.
///
/// ## Example
///
/// ```
/// # use gsl_splines::{Accelerator, Interpolation2dType, Spline2d, SplineError};
/// # use ndarray::{Array1, Array2};
/// #
/// # fn main() -> Result<(), SplineError> {
/// // Data creation
/// let xdata = Array1::linspace(0.0, 5.0 ,10);
/// let ydata = Array1::linspace(0.0, 8.0 ,20);
/// let zdata = ndarray::Array2::<f64>::ones((xdata.len(), ydata.len()));
///
/// // Interpolation type and Accelerator
/// let xacc = gsl_splines::Accelerator::new();
/// let yacc = gsl_splines::Accelerator::new();
/// let typ = gsl_splines::Interpolation2dType::Bicubic;
///
/// // Spline creation
/// let spline2d = gsl_splines::Spline2d::build(typ, xdata, ydata, zdata, xacc, yacc);
///
/// println!("{:#?}", spline2d);
/// # Ok(())
/// # }
/// ```
/// [`type`]: enum.InterpolationType.html
/// [`Accelerators`]: struct.Accelerator.html
#[allow(dead_code)]
pub struct Spline2d {
    /// Interpolation Type.
    pub typ: Interpolation2dType,
    /// Array of x data points.
    pub xdata: Array1<f64>,
    /// Array of y data points.
    pub ydata: Array1<f64>,
    /// Grid of z data points.
    // TODO: Consider if storing the 2d array is necessary, since into_flatten() could skip the copying.
    pub zdata: Array2<f64>,
    /// Size of x array.
    pub(crate) zdata_flat: Array1<f64>,
    pub(crate) xsize: usize,
    /// Size of y array.
    pub(crate) ysize: usize,
    /// `rgsl`s Spline object
    pub(crate) shape: (usize, usize),
    pub(crate) gsl_spline2d: RgslSpline2d,
    /// Reference to a common [`Accelerator`] for the x data points.
    ///
    /// [`Accelerator`]: struct.Accelerator.html
    pub xacc: Rc<RefCell<Accelerator>>,
    /// Reference to a common [`Accelerator`] for the y data points.
    ///
    /// [`Accelerator`]: struct.Accelerator.html
    pub yacc: Rc<RefCell<Accelerator>>,
}

impl Spline2d {
    /// Creates a new Spline2d with common `xacc` and `yacc` [`Accelerators`].
    ///
    /// [`Accelerators`]: struct.Accelerator.html
    pub fn build(
        typ: Interpolation2dType,
        xdata: Array1<f64>,
        ydata: Array1<f64>,
        zdata: Array2<f64>,
        xacc: Rc<RefCell<Accelerator>>,
        yacc: Rc<RefCell<Accelerator>>,
    ) -> Result<Self> {
        Spline2d::check_data(&xdata, &ydata, &zdata, typ)?;
        let xdata = xdata.clone();
        let ydata = ydata.clone();
        let zdata = zdata.clone();
        let zdata_flat = zdata
            .flatten_with_order(ndarray::Order::RowMajor)
            .to_owned();
        let xsize = xdata.len();
        let ysize = ydata.len();
        let shape = (xsize, ysize);

        // Alloc the gsl_spline2d, initialize it later
        let gsl_spline2d = Spline2d::build_gsl_spline2d(typ.into(), xsize, ysize)?;

        let mut spline2d = Spline2d {
            typ,
            xdata,
            ydata,
            zdata,
            zdata_flat,
            xsize,
            ysize,
            shape,
            gsl_spline2d,
            xacc,
            yacc,
        };
        spline2d.init_gsl_spline2d()?;
        Ok(spline2d)
    }

    /// Checks if supplied datasets are valid.
    fn check_data(
        x: &Array1<f64>,
        y: &Array1<f64>,
        z: &Array2<f64>,
        typ: Interpolation2dType,
    ) -> Result<()> {
        if x.len() <= 1 || y.len() <= 1 {
            return Err(SplineError::InvalidDataset);
        }
        let gsl_typ: RgslInterp2dType = typ.into();
        if x.len() + y.len() < gsl_typ.min_size() as usize {
            return Err(SplineError::NotEnoughPoints);
        }
        if !x.iter().is_sorted() {
            return Err(SplineError::UnsortedDataset("x".into()));
        }
        if !y.iter().is_sorted() {
            return Err(SplineError::UnsortedDataset("y".into()));
        }

        if x.len() * y.len() != z.len() {
            return Err(SplineError::Dataset2dMismatch);
        }
        if (x.len() != z.shape()[0]) | (y.len() != z.shape()[1]) {
            return Err(SplineError::Dataset2dMismatch);
        }
        Ok(())
    }

    /// Creates a new uninitialized `rgsl::Spline2d` object and returns it.
    fn build_gsl_spline2d(
        typ: RgslInterp2dType,
        xsize: usize,
        ysize: usize,
    ) -> Result<RgslSpline2d> {
        // Calls gsl_spline2d_alloc()
        RgslSpline2d::new(typ, xsize, ysize).ok_or(SplineError::GSLInterp2dAlloc)
    }

    /// Initializes the interpolation object
    fn init_gsl_spline2d(&mut self) -> Result<()> {
        // Calls gsl_spline_init()
        // Safe unwrap, arrays are already checked
        let temp = self.gsl_spline2d.init(
            self.xdata.view().as_standard_layout().as_slice().unwrap(),
            self.ydata.view().as_standard_layout().as_slice().unwrap(),
            self.zdata.view().as_standard_layout().as_slice().unwrap(),
        );
        match temp {
            Ok(_) => Ok(()),
            Err(err) => Err(SplineError::GSLSpline2dInit { err }),
        }
    }

    /// Retruns the [`Interpolation2dType`]'s name.
    ///
    /// [`Interpolation2dType`]: enum.InterpolationType.html
    pub fn name(&self) -> String {
        self.gsl_spline2d.name()
    }

    /// Returns the [`Interpolation2dType`]'s minimum required number of x points. For example,
    /// bicubic interpolation requires a minimum of 4 points.
    ///
    /// [`Interpolation2dType`]: enum.InterpolationType.html
    pub fn min_size(&self) -> usize {
        self.gsl_spline2d.min_size()
    }

    /// Resets both [`Accelerator`]s' cache and stats.
    ///
    /// [`Accelerator`]: struct.Accelerator.html
    pub fn reset_acc(&mut self) {
        self.yacc.borrow_mut().reset();
        self.xacc.borrow_mut().reset();
    }

    /// Sets the value z_ij of grid point (i, j) of the stored zdata array to z.
    pub fn set(&mut self, i: usize, j: usize, z: f64) {
        self.gsl_spline2d
            .set(self.zdata_flat.as_slice_mut().unwrap(), i, j, z);
        self.zdata[[i, j]] = z;
    }

    /// Returns the value z_ij of grid point (i, j) of the stored zdata array
    pub fn get(&mut self, i: usize, j: usize) -> f64 {
        self.gsl_spline2d
            .get(self.zdata_flat.as_slice_mut().unwrap(), i, j)
    }

    /// Returns the index corresponding to the grid point (i, j), which is given by
    /// `j*xsize + i`.
    pub fn idx(&mut self, i: usize, j: usize) -> usize {
        j * self.xsize + i
    }
}

impl std::fmt::Debug for Spline2d {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Spline")
            .field("Interpolation Type", &self.typ)
            .field(
                "xdata",
                &format_args!(
                    "[{:.5}, .. {:.5}], size={}",
                    self.xdata[0],
                    self.xdata[self.xsize - 1],
                    self.xsize
                ),
            )
            .field(
                "ydata",
                &format_args!(
                    "[{:.5}, .. {:.5}], size={}",
                    self.ydata[0],
                    self.ydata[self.ysize - 1],
                    self.ysize
                ),
            )
            .field("zdata", &self.zdata.to_string())
            .field("X-Accelerator", &self.xacc.borrow_mut())
            .field("X-Accelerator", &self.yacc.borrow_mut())
            .finish()
        // writeln!(f, "\nzdata:\n{}", self.zdata)
    }
}
