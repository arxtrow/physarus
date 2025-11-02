use std::ops::*;

use opimps::impl_ops;

use crate::core::{
    traits::*,
    vec_states::{ScaleSeal, SpaceSeal},
    *,
};

pub type Mat2 = Matrix<2, 2>;

impl Mat2 {
    #[rustfmt::skip]
    pub fn new([
        m00, m01,
        m10, m11,
    ]: [Float; 4]) -> Self {
        Self {
            cols: [
                [m00, m10],
                [m01, m11],
            ]
        }
    }

    pub fn det(&self) -> Float {
        let [[m00, m10], [m01, m11]] = self.to_cols();

        m00 * m11 - m01 * m10
    }

    pub fn inverse(&self) -> Option<Self> {
        let [[m00, m10], [m01, m11]] = self.to_cols();

        let det = self.det();

        if det.is_zero_eps() {
            return None;
        }

        let inv_det = 1.0 / det;

        let c00 = m11;
        let c01 = -m10;
        let c10 = -m01;
        let c11 = m00;

        let inverse = inv_det * Self::from_cols([[c00, c01], [c10, c11]]);

        Some(inverse)
    }

    pub fn from_angle(angle: Float) -> Self {
        let (sin, cos) = angle.sin_cos();

        Self::from_cols([[cos, sin], [-sin, cos]])
    }

    pub fn from_scale_array(scale: [Float; 2]) -> Self {
        Self::from_cols([[scale[0], 0.0], [0.0, scale[1]]])
    }

    pub fn from_scale(scale: Float) -> Self {
        Self::from_cols([[scale, 0.0], [0.0, scale]])
    }
}

// Specialized Mat2 * Mat2
#[impl_ops(Mul)]
fn mul(self: Mat2, rhs: Mat2) -> Mat2 {
    let [[a00, a10], [a01, a11]] = self.to_cols();
    let [[b00, b10], [b01, b11]] = rhs.to_cols();

    Mat2::from_cols([
        [a00 * b00 + a01 * b10, a10 * b00 + a11 * b10],
        [a00 * b01 + a01 * b11, a10 * b01 + a11 * b11],
    ])
}

// Specialized Mat2 * Vector
#[impl_ops(Mul)]
fn mul<Sem, Scale: ScaleSeal, Space: SpaceSeal>(
    self: Mat2,
    vec: Vector<2, Sem, Scale, Space>,
) -> Direction2 {
    Direction::from_array_in([
        self[(0, 0)] * vec[0] + self[(0, 1)] * vec[1],
        self[(1, 0)] * vec[0] + self[(1, 1)] * vec[1],
    ])
}

/// Creates a 2x2 matrix with custom formatting.
#[macro_export]
macro_rules! mat2 {
    {$($values:expr),* $(,)?} => {
        Mat2::new([$($values),*])
    };
}

pub use mat2;
