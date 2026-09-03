//! Definitions of `Machine` and `MachineBuilder`.

use crate::{Bfield, Current, Geometry, Perturbation, Qfactor};

/// To be used when no perturbation is passed, to avoid calling `unwrap()` every time.
const ZERO_PERTURBATION: &Perturbation = &Perturbation::zero();

/// Used to create a [`Machine`].
///
/// This object only holds **references** to the machine objects.
#[derive(Debug, Clone, Copy)]
#[expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]
pub struct MachineBuilder<'obj> {
    geometry: Option<&'obj dyn Geometry>,
    qfactor: &'obj dyn Qfactor,
    current: &'obj dyn Current,
    bfield: &'obj dyn Bfield,
    perturbation: Option<&'obj Perturbation>,
}

impl<'obj> MachineBuilder<'obj> {
    /// Creates a new `MachineBuilder` from the machine's q-factor, plasma current and magnetic field.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    ///
    /// let builder = MachineBuilder::new(&qfactor, &current, &bfield);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn new(
        qfactor: &'obj dyn Qfactor,
        current: &'obj dyn Current,
        bfield: &'obj dyn Bfield,
    ) -> Self {
        Self {
            geometry: None,
            qfactor,
            current,
            bfield,
            perturbation: None,
        }
    }

    /// Sets the machine's [`Geometry`] object.
    ///
    /// The [`Geometry`] object is only needed for conversions to and from laboratory coordinates
    /// and is not used in any routines.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let geometry = LarGeometry::new(2.5, 1.75, 0.5);
    /// let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    ///
    /// let builder = MachineBuilder::new(&qfactor, &current, &bfield).with_geometry(&geometry);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn with_geometry(mut self, geometry: &'obj dyn Geometry) -> Self {
        self.geometry = Some(geometry);
        self
    }

    /// Sets the machine's [`Perturbation`] object.
    ///
    /// If no perturbation is added, then it is set to [`Perturbation::zero()`].
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let qfactor = UnityQfactor::new(lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    ///
    /// let builder = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_perturbation(&perturbation);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn with_perturbation(mut self, perturbation: &'obj Perturbation) -> Self {
        self.perturbation = Some(perturbation);
        self
    }

    /// Creates a new [`Machine`] with the Builder's configuration.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let geometry = LarGeometry::new(2.5, 1.75, 0.5);
    /// let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    ///
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_geometry(&geometry)
    ///     .with_perturbation(&perturbation)
    ///     .build();
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn build(self) -> Machine<'obj> {
        Machine::build(self)
    }
}

/// Container for all information about the magnetic field, geometry and perturbations of a device.
///
/// This object only holds references to the machine objects.
///
/// Should be created with a [`MachineBuilder`].
#[derive(Debug, Clone, Copy)]
#[expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]
pub struct Machine<'obj> {
    geometry: Option<&'obj dyn Geometry>,
    qfactor: &'obj dyn Qfactor,
    current: &'obj dyn Current,
    bfield: &'obj dyn Bfield,
    perturbation: Option<&'obj Perturbation>,
}

impl<'obj> Machine<'obj> {
    /// Constructs a `Machine` from a [`MachineBuilder`].
    #[must_use]
    pub(crate) fn build(builder: MachineBuilder<'obj>) -> Self {
        Machine {
            geometry: builder.geometry,
            qfactor: builder.qfactor,
            current: builder.current,
            bfield: builder.bfield,
            perturbation: builder.perturbation,
        }
    }

    /// Returns a reference to the [`Geometry`] object.
    #[must_use]
    pub fn geometry(&self) -> Option<&dyn Geometry> {
        self.geometry
    }

    /// Returns a reference to the [`Qfactor`] object.
    #[must_use]
    pub fn qfactor(&self) -> &dyn Qfactor {
        self.qfactor
    }

    /// Returns a reference to the [`Current`] object.
    #[must_use]
    pub fn current(&self) -> &dyn Current {
        self.current
    }

    /// Returns a reference to the [`Bfield`] object.
    #[must_use]
    pub fn bfield(&self) -> &dyn Bfield {
        self.bfield
    }

    /// Returns a reference to the [`Perturbation`] object.
    #[must_use]
    pub fn perturbation(&self) -> &Perturbation {
        // To avoid calling `unwrap` every time
        self.perturbation.unwrap_or(ZERO_PERTURBATION)
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn machine_init() {
        let geometry = LarGeometry::new(2.5, 1.75, 0.5);
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let perturbation = Perturbation::zero();

        let machine = MachineBuilder::new(&qfactor, &current, &bfield)
            .with_geometry(&geometry)
            .with_perturbation(&perturbation)
            .build();

        let _: Option<&dyn Geometry> = machine.geometry();
        let _: &dyn Qfactor = machine.qfactor();
        let _: &dyn Current = machine.current();
        let _: &dyn Bfield = machine.bfield();
        let _: &Perturbation = machine.perturbation();
    }
}
