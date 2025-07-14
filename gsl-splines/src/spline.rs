use crate::{Accelerator, InterpolationType, Result, RgslInterpType, RgslSpline, SplineError};

use ndarray::Array1;
use std::cell::RefCell;
use std::rc::Rc;

/// A Spline of a specific type that stores the data points and can share an Accelerator with other
/// splines.
///
/// ## Note
///
/// `gsl_interp_free()` is called inside `rgsl`.
pub struct Spline {
    /// Interpolation Type.
    pub typ: InterpolationType,
    /// Array of x data points.
    pub xdata: Array1<f64>,
    /// Array of y data points.
    pub ydata: Array1<f64>,
    /// Size of x and y arrays.
    pub size: usize,
    /// `rgsl`s Spline object
    pub(crate) gsl_spline: RgslSpline,
    /// Reference to a common `Accelerator` that can be used by many splines that are to be evaluated
    /// at the same point.
    pub acc: Rc<RefCell<Accelerator>>,
}

impl Spline {
    /// Creates a new `Spline` with a common `Accelerator`
    pub fn build(
        typ: InterpolationType,
        xdata: &Array1<f64>,
        ydata: &Array1<f64>,
        acc: Rc<RefCell<Accelerator>>,
    ) -> Result<Self> {
        Spline::check_data(xdata, ydata, typ)?;
        let xdata = xdata.clone();
        let ydata = ydata.clone();
        let size = xdata.len();

        // Alloc the gsl_spline, initialize it later
        let gsl_spline = Spline::build_gsl_spline(typ.into(), size)?;

        let mut interp = Spline {
            typ,
            xdata,
            ydata,
            size,
            gsl_spline,
            acc,
        };

        interp.init_gsl_interp()?;
        Ok(interp)
    }

    /// Checks if supplied datasets are valid.
    fn check_data(x: &Array1<f64>, y: &Array1<f64>, typ: InterpolationType) -> Result<()> {
        if x.len() <= 1 || y.len() <= 1 {
            return Err(SplineError::InvalidDataset);
        }
        let gsl_typ: RgslInterpType = typ.into();
        if x.len() < gsl_typ.min_size() as usize {
            return Err(SplineError::NotEnoughPoints);
        }
        if x.len() != y.len() {
            return Err(SplineError::DatasetMismatch);
        }
        if !x.iter().is_sorted() {
            return Err(SplineError::UnsortedDataset);
        }
        Ok(())
    }

    /// Creates a new uninitialized `rgsl::Interp` object and returns it.
    fn build_gsl_spline(typ: RgslInterpType, size: usize) -> Result<RgslSpline> {
        // Calls gsl_spline_alloc()
        RgslSpline::new(typ, size).ok_or(SplineError::GSLInterpAlloc)
    }

    /// Initializes the interpolation object
    fn init_gsl_interp(&mut self) -> Result<()> {
        // Calls gsl_spline_init()
        // Safe unwrap, arrays are already checked
        let temp = self.gsl_spline.init(
            self.xdata.view().as_standard_layout().as_slice().unwrap(),
            self.ydata.view().as_standard_layout().as_slice().unwrap(),
        );
        match temp {
            Ok(_) => Ok(()),
            Err(err) => Err(SplineError::GSLSplineInit { err }),
        }
    }

    /// Retruns the Interpolation type's name.
    pub fn name(&self) -> String {
        self.gsl_spline.name()
    }

    /// Returns the interpolation type's minimum required number of x points.
    pub fn min_size(&self) -> usize {
        self.gsl_spline.min_size() as usize
    }

    /// Resets the `Accelerator`'s cache and stats.
    pub fn reset_acc(&mut self) {
        self.acc.borrow_mut().reset();
    }
}

impl std::fmt::Debug for Spline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Spline")
            .field("Interpolation Type", &self.typ)
            .field(
                "xdata",
                &format_args!(
                    "[{:.5}, .. {:.5}], size={}",
                    self.xdata[0],
                    self.xdata[self.size - 1],
                    self.size
                ),
            )
            .field(
                "ydata",
                &format_args!(
                    "[{:.5}, .. {:.5}], size={}",
                    self.ydata[0],
                    self.ydata[self.size - 1],
                    self.size
                ),
            )
            .field("Accelerator", &self.acc.borrow_mut())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn some_spline() -> Spline {
        let xdata = Array1::linspace(0.0, 3.0, 10);
        let ydata = xdata.pow2();
        let typ = InterpolationType::Akima;
        let acc = Accelerator::new();
        Spline::build(typ, &xdata, &ydata, acc).unwrap()
    }

    #[test]
    fn test_check_data() {
        // Also check that error messages print
        assert!(matches!(
            Spline::check_data(
                &Array1::from_vec(vec![1.0]),
                &Array1::from_vec(vec![2.0]),
                InterpolationType::Linear,
            )
            .unwrap_err(),
            SplineError::InvalidDataset
        ));
        let _ = format!("{:?}", SplineError::InvalidDataset);

        assert!(matches!(
            Spline::check_data(
                &Array1::from_vec(vec![1.0, 2.0]),
                &Array1::from_vec(vec![2.0, 3.0]),
                InterpolationType::Akima,
            )
            .unwrap_err(),
            SplineError::NotEnoughPoints
        ));
        let _ = format!("{:?}", SplineError::NotEnoughPoints);

        assert!(matches!(
            Spline::check_data(
                &Array1::from_vec(vec![1.0, 2.0]),
                &Array1::from_vec(vec![2.0, 3.0, 4.0]),
                InterpolationType::Linear,
            )
            .unwrap_err(),
            SplineError::DatasetMismatch
        ));
        let _ = format!("{:?}", SplineError::DatasetMismatch);

        assert!(matches!(
            Spline::check_data(
                &Array1::from_vec(vec![2.0, 1.0]),
                &Array1::from_vec(vec![2.0, 3.0]),
                InterpolationType::Linear,
            )
            .unwrap_err(),
            SplineError::UnsortedDataset
        ));
        let _ = format!("{:?}", SplineError::UnsortedDataset);
    }

    #[test]
    fn test_spline_build() {
        let xdata = Array1::linspace(0.0, 3.0, 10);
        let ydata = xdata.pow2();
        let typ = InterpolationType::Akima;
        let acc = Accelerator::new();
        assert!(Spline::build(typ, &xdata, &ydata, acc).is_ok());
    }

    #[test]
    fn test_name() {
        let spline = some_spline();
        spline.name();
    }

    #[test]
    fn test_min_size() {
        let spline = some_spline();
        spline.min_size();
    }

    #[test]
    fn test_reset_acc() {
        let mut spline = some_spline();
        spline.reset_acc();
    }

    #[test]
    fn test_debug_trait() {
        let spline = some_spline();
        let _ = format!("{spline:?}");
        let _ = format!("{spline:#?}");
    }
}
