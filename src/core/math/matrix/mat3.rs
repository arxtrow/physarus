use std::{array, ops::*, simd::Simd};

use opimps::impl_ops;

use crate::core::{traits::*, vec_states::*, *};

pub type Mat3 = Matrix<3, 3>;

impl Mat3 {
    #[rustfmt::skip]
    pub fn new([
        m00, m01, m02,
        m10, m11, m12,
        m20, m21, m22,
    ]: [Float; 9]) -> Self {
        Self {
            cols: [
                [m00, m10, m20],
                [m01, m11, m21],
                [m02, m12, m22],
            ]
        }
    }

    pub fn det(&self) -> Float {
        let [[m00, m10, m20], [m01, m11, m21], [m02, m12, m22]] = self.to_cols();

        let u0 = m00;
        let l12 = (m11, m22).diff_of_mul((m12, m21));

        let u1 = m01;
        let l02 = (m10, m22).diff_of_mul((m12, m20));

        let u2 = m02;
        let l01 = (m10, m21).diff_of_mul((m11, m20));

        // Sum over C(3,1)=3 ways to partition columns into 1D×2D orthogonal subspaces.
        // Signs: [0,1,2]→+, [1,0,2]→-, [2,0,1]→+ (counting index swaps)
        u0 * l12 - u1 * l02 + u2 * l01
    }

    pub fn inverse(&self) -> Option<Self> {
        let [[m00, m10, m20], [m01, m11, m21], [m02, m12, m22]] = self.to_cols();

        let u0 = m00;
        let l12 = (m11, m22).diff_of_mul((m12, m21));

        let u1 = m01;
        let l02 = (m10, m22).diff_of_mul((m12, m20));

        let u2 = m02;
        let l01 = (m10, m21).diff_of_mul((m11, m20));

        let det = u0 * l12 - u1 * l02 + u2 * l01;

        if det.is_zero_eps() {
            return None;
        }

        let inv_det = 1.0 / det;

        let c00 = l12;
        let c01 = -l02;
        let c02 = l01;

        let c10 = -(m01, m22).diff_of_mul((m02, m21));
        let c11 = (m00, m22).diff_of_mul((m02, m20));
        let c12 = -(m00, m21).diff_of_mul((m01, m20));

        let c20 = (m01, m12).diff_of_mul((m02, m11));
        let c21 = -(m00, m12).diff_of_mul((m02, m10));
        let c22 = (m00, m11).diff_of_mul((m01, m10));

        let inverse =
            inv_det * Self::from_cols([[c00, c01, c02], [c10, c11, c12], [c20, c21, c22]]);

        Some(inverse)
    }

    pub fn to_hom_4x4(&self) -> Mat4 {
        let [[m00, m10, m20], [m01, m11, m21], [m02, m12, m22]] = self.to_cols();

        mat4!(
            m00, m01, m02, 0.0,
            m10, m11, m12, 0.0,
            m20, m21, m22, 0.0,
            0.0, 0.0, 0.0, 1.0,
        )
    }
}

// Specialized Mat3 * Mat3
#[impl_ops(Mul)]
fn mul(self: Mat3, rhs: Mat3) -> Mat3 {
    let a = self;

    let a_col0 = a.col(0);
    let a_col1 = a.col(1);
    let a_col2 = a.col(2);

    let a_col0 = Simd::from_array([a_col0[0], a_col0[1], a_col0[2], 0.0]);
    let a_col1 = Simd::from_array([a_col1[0], a_col1[1], a_col1[2], 0.0]);
    let a_col2 = Simd::from_array([a_col2[0], a_col2[1], a_col2[2], 0.0]);

    let mut c_cols = [[0.0; 3]; 3];

    for j in 0..3 {
        let b = rhs.col(j);
        let b0j = Simd::splat(b[0]);
        let b1j = Simd::splat(b[1]);
        let b2j = Simd::splat(b[2]);

        let c_col_j = [b0j, b1j, b2j].linear_combination([a_col0, a_col1, a_col2]);

        let res_array = c_col_j.to_array();

        c_cols[j] = [res_array[0], res_array[1], res_array[2]];
    }

    Mat3::from_cols(c_cols)
}

// Specialized Mat3 * Vector
#[impl_ops(Mul)]
fn mul<Sem, Scale: ScaleSeal, Space: SpaceSeal>(
    self: Mat3,
    vec: Vec3<Sem, Scale, Space>,
) -> Direction3 {
    let a = self;

    let a_col0 = a.col(0);
    let a_col1 = a.col(1);
    let a_col2 = a.col(2);

    let a_col0 = Simd::from_array([a_col0[0], a_col0[1], a_col0[2], 0.0]);
    let a_col1 = Simd::from_array([a_col1[0], a_col1[1], a_col1[2], 0.0]);
    let a_col2 = Simd::from_array([a_col2[0], a_col2[1], a_col2[2], 0.0]);

    let v = vec.as_slice();
    let v0 = Simd::splat(v[0]);
    let v1 = Simd::splat(v[1]);
    let v2 = Simd::splat(v[2]);

    let w = [v0, v1, v2].linear_combination([a_col0, a_col1, a_col2]);
    let [x, y, z, _] = w.to_array();

    Direction3::new(x, y, z)
}

/// Creates a 3x3 matrix with custom formatting.
#[macro_export]
macro_rules! mat3 {
    {$($values:expr),* $(,)?} => {
        Mat3::new([$($values),*])
    };
}

pub use mat3;
