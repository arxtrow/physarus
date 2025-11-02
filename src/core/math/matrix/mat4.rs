use std::{ops::*, simd::Simd};

use opimps::impl_ops;

use crate::prelude::*;

pub type Mat4 = Matrix<4, 4>;

impl Mat4 {
    #[rustfmt::skip]
    pub fn new([
        m00, m01, m02, m03,
        m10, m11, m12, m13,
        m20, m21, m22, m23,
        m30, m31, m32, m33,
    ]: [Float; 16]) -> Self {
        Self {
            cols: [
                [m00, m10, m20, m30],
                [m01, m11, m21, m31],
                [m02, m12, m22, m32],
                [m03, m13, m23, m33],
            ]
        }
    }

    #[inline]
    pub fn upper3x3(&self) -> Mat3 {
        self.sub_mat(0..3, 0..3)
    }

    pub fn det(&self) -> Float {
        let [
            [m00, m10, m20, m30],
            [m01, m11, m21, m31],
            [m02, m12, m22, m32],
            [m03, m13, m23, m33],
        ] = self.to_cols();

        //
        // Upper(u) block (rows 0,1):  [m00  m01  m02  m03]
        //                             [m10  m11  m12  m13]
        //
        // Lower(l) block (rows 2,3):  [m20  m21  m22  m23]
        //                             [m30  m31  m32  m33]
        //
        //

        // Bivector areas (signed 2D projections) for each column pair.
        // u_ij = area of columns {i,j} projected onto upper plane (e0∧e1)
        // l_ij = area of columns {i,j} projected onto lower plane (e2∧e3)
        let u01 = (m00, m11).diff_of_mul((m01, m10));
        let l23 = (m22, m33).diff_of_mul((m23, m32));

        let u02 = (m00, m12).diff_of_mul((m02, m10));
        let l13 = (m21, m33).diff_of_mul((m23, m31));

        let u03 = (m00, m13).diff_of_mul((m03, m10));
        let l12 = (m21, m32).diff_of_mul((m22, m31));

        let u12 = (m01, m12).diff_of_mul((m02, m11));
        let l03 = (m20, m33).diff_of_mul((m23, m30));

        let u13 = (m01, m13).diff_of_mul((m03, m11));
        let l02 = (m20, m32).diff_of_mul((m22, m30));

        let u23 = (m02, m13).diff_of_mul((m03, m12));
        let l01 = (m20, m31).diff_of_mul((m21, m30));

        // Each u_ij × l_kl gives 4D volume of |u_ij|(e0^e1) ^ |l_kl|(e2^e3).
        //
        // For example:
        // u01×l23 measures the 4D volume when:
        //   - Columns {0,1} span a parallelogram in the e0^e1 plane with signed area u01
        //   - Columns {2,3} span a parallelogram in the e2^e3 plane with signed area l23
        //
        // u03×l12 measures the 4D volume when:
        //   - Columns {0,3} span a parallelogram in the e0^e1 plane with signed area u03
        //   - Columns {1,2} span a parallelogram in the e2^e3 plane with signed area l12
        //
        // We sum over all C(4,2) = 6 ways to partition 4 columns into two disjoint pairs.
        // Each partition splits the 4D space differently into orthogonal 2D×2D subspaces.
        //
        // Signs align to canonical basis e0^e1^e2^e3
        //   [0,1,2,3]→(+)  [0,2,1,3]→(-)  [0,3,1,2]→(+)
        //   [1,2,0,3]→(+)  [1,3,0,2]→(-)  [2,3,0,1]→(+)
        u01 * l23 - u02 * l13 + u03 * l12 + u12 * l03 - u13 * l02 + u23 * l01
    }

    pub fn inverse(&self) -> Option<Self> {
        let [
            [m00, m10, m20, m30],
            [m01, m11, m21, m31],
            [m02, m12, m22, m32],
            [m03, m13, m23, m33],
        ] = self.to_cols();

        let u01 = (m00, m11).diff_of_mul((m01, m10));
        let l23 = (m22, m33).diff_of_mul((m23, m32));

        let u02 = (m00, m12).diff_of_mul((m02, m10));
        let l13 = (m21, m33).diff_of_mul((m23, m31));

        let u03 = (m00, m13).diff_of_mul((m03, m10));
        let l12 = (m21, m32).diff_of_mul((m22, m31));

        let u12 = (m01, m12).diff_of_mul((m02, m11));
        let l03 = (m20, m33).diff_of_mul((m23, m30));

        let u13 = (m01, m13).diff_of_mul((m03, m11));
        let l02 = (m20, m32).diff_of_mul((m22, m30));

        let u23 = (m02, m13).diff_of_mul((m03, m12));
        let l01 = (m20, m31).diff_of_mul((m21, m30));

        let det = u01 * l23 - u02 * l13 + u03 * l12 + u12 * l03 - u13 * l02 + u23 * l01;

        if det.is_zero_eps() {
            return None;
        }

        let inv_det = 1.0 / det;

        let c00 = (m11 * l23 - m12 * l13 + m13 * l12);
        let c01 = -(m10 * l23 - m12 * l03 + m13 * l02);
        let c02 = (m10 * l13 - m11 * l03 + m13 * l01);
        let c03 = -(m10 * l12 - m11 * l02 + m12 * l01);

        let c10 = -(m01 * l23 - m02 * l13 + m03 * l12);
        let c11 = (m00 * l23 - m02 * l03 + m03 * l02);
        let c12 = -(m00 * l13 - m01 * l03 + m03 * l01);
        let c13 = (m00 * l12 - m01 * l02 + m02 * l01);

        let c20 = (m31 * u23 - m32 * u13 + m33 * u12);
        let c21 = -(m30 * u23 - m32 * u03 + m33 * u02);
        let c22 = (m30 * u13 - m31 * u03 + m33 * u01);
        let c23 = -(m30 * u12 - m31 * u02 + m32 * u01);

        let c30 = -(m21 * u23 - m22 * u13 + m23 * u12);
        let c31 = (m20 * u23 - m22 * u03 + m23 * u02);
        let c32 = -(m20 * u13 - m21 * u03 + m23 * u01);
        let c33 = (m20 * u12 - m21 * u02 + m22 * u01);

        let cof_t = mat4!( // transposed
            c00, c10, c20, c30,
            c01, c11, c21, c31,
            c02, c12, c22, c32,
            c03, c13, c23, c33,
        );

        let inverse = inv_det * cof_t;

        Some(inverse)
    }
}

#[impl_ops(Mul)]
fn mul(self: Mat4, rhs: Mat4) -> Mat4 {
    let a = self;
    let b = rhs;

    let a_col0 = Simd::from_array(a.col(0));
    let a_col1 = Simd::from_array(a.col(1));
    let a_col2 = Simd::from_array(a.col(2));
    let a_col3 = Simd::from_array(a.col(3));

    let mut c_cols = [[0.0; 4]; 4];

    for j in 0..4 {
        let b_col = b.col(j);

        let b0j = Simd::splat(b_col[0]);
        let b1j = Simd::splat(b_col[1]);
        let b2j = Simd::splat(b_col[2]);
        let b3j = Simd::splat(b_col[3]);

        // For C = A * B
        // Column j of C is a linear combination of A's columns, weighted by coeffs of column j of B:
        // C[:, j] = B[0, j]*A[:, 0] + B[1, j]*A[:, 1] + B[2, j]*A[:, 2] + B[3, j]*A[:, 3]
        let c_col_j = [b0j, b1j, b2j, b3j].linear_combination([a_col0, a_col1, a_col2, a_col3]);

        c_cols[j] = c_col_j.to_array();
    }

    Mat4::from_cols(c_cols)
}

#[impl_ops(Mul)]
fn mul<Sem, Scale: ScaleSeal, Space: SpaceSeal>(
    self: Mat4,
    vec: Vec4<Sem, Scale, Space>,
) -> Direction4 {
    let a = self;

    let a_col0 = Simd::from_array(a.col(0));
    let a_col1 = Simd::from_array(a.col(1));
    let a_col2 = Simd::from_array(a.col(2));
    let a_col3 = Simd::from_array(a.col(3));

    let v = vec.as_slice();
    let v0 = Simd::splat(v[0]);
    let v1 = Simd::splat(v[1]);
    let v2 = Simd::splat(v[2]);
    let v3 = Simd::splat(v[3]);

    // w = A * v is a linear combination of A's columns, weighted by components of v:
    // w = v[0]*A[:,0] + v[1]*A[:,1] + v[2]*A[:,2] + v[3]*A[:,3]
    let w = [v0, v1, v2, v3].linear_combination([a_col0, a_col1, a_col2, a_col3]);

    Direction::from_array_in(w.to_array())
}

/// Creates a 4x4 matrix with custom formatting.
///
/// This macro disables rustfmt rules, allowing you to manually align matrix entries
/// for better readability. Without it, `Mat4::new()` would format all entries on one line.
/// # Examples
///
/// ```
/// use physarus::{mat4, prelude::Mat4};
///
/// let long_ass_name = 1.0;
///
/// let m = mat4!(
///     1.0,   2.0,                 323.0, 4.0,
///     2.234, 3.0,                 8.0,   7.0,
///     5.0,   4.0 / long_ass_name, 2.0,   1.0,
///     0.0,   0.0,                 0.0,   1.0,
/// );
///
/// ```
#[macro_export]
macro_rules! mat4 {
    {$($values:expr),* $(,)?} => {
        Mat4::new([$($values),*])
    };
}

pub use mat4;
