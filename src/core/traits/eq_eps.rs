use crate::prelude::*;

#[cfg(not(feature = "f64"))]
const EPS: Float = 1e-6; // Float::EPSILON too small

#[cfg(feature = "f64")]
const EPS: Float = 1e-12;

pub trait EqEps {
    /// Compares two floats for approximate equality.
    ///
    /// # Returns
    /// - `true` if the values are "close enough"; `false` otherwise.
    fn eq_eps(&self, b: &Self) -> bool;

    /// Checks if a float is "close enough to zero" as defined by the current EPS.
    fn is_zero_eps(&self) -> bool;
}

impl EqEps for Float {
    fn eq_eps(&self, lhs: &Self) -> bool {
        let rhs = *self;
        let lhs = *lhs;

        match (rhs, lhs) {
            (rhs, lhs) if rhs.is_nan() || lhs.is_nan() => false,

            (rhs, lhs) if rhs == lhs => true,

            (rhs, lhs) if rhs.abs().is_zero_eps() && lhs.abs().is_zero_eps() => true,

            (rhs, lhs) => {
                let diff = (rhs - lhs).abs();

                let relative_epsilon = EPS * (rhs.abs() + lhs.abs()) * 0.5;
                diff <= relative_epsilon.max(EPS)
            }
        }
    }

    fn is_zero_eps(&self) -> bool {
        self.abs() <= EPS
    }
}
