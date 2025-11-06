use crate::prelude::*;
use std::{array, ops::Range};

/// Matrix stored in column-major order: array of column vectors
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Matrix<const R: usize, const C: usize> {
    pub(in crate::core) cols: [[Float; R]; C],
}

impl<const R: usize, const C: usize> Matrix<R, C> {
    /// Create matrix from column vectors. ***Prefered way***
    #[inline]
    pub fn from_cols(cols: [[Float; R]; C]) -> Self {
        Self { cols }
    }

    /// Create matrix from row vectors
    #[inline]
    pub fn from_rows(rows: [[Float; C]; R]) -> Self {
        let cols = array::from_fn(|col| array::from_fn(|row| rows[row][col]));

        Self { cols }
    }

    /// Create matrix from column-major array (flat array in column-major layout)
    #[inline]
    pub fn from_cols_flat(data: &[Float]) -> Self {
        debug_assert_eq!(data.len(), R * C);

        let cols = array::from_fn(|col| array::from_fn(|row| data[col * R + row]));

        Self { cols }
    }

    /// Create matrix from row-major array (flat array in row-major layout)
    #[inline]
    pub fn from_rows_flat(data: &[Float]) -> Self {
        debug_assert_eq!(data.len(), R * C);

        let cols = array::from_fn(|col| array::from_fn(|row| data[row * C + col]));

        Self { cols }
    }

    /// Create zero matrix
    #[inline]
    pub const fn zero() -> Self {
        Self {
            cols: [[0.0; R]; C],
        }
    }

    /// Get column as new vector
    #[inline]
    pub fn col(&self, col: usize) -> [Float; R] {
        debug_assert!(col < C);
        self.cols[col]
    }

    /// Get mutable reference to column
    #[inline]
    pub fn col_mut(&mut self, col: usize) -> &mut [Float; R] {
        debug_assert!(col < C);
        &mut self.cols[col]
    }

    /// Get row as new vector
    #[inline]
    pub fn row(&self, row: usize) -> [Float; C] {
        debug_assert!(row < R);
        array::from_fn(|col| self.cols[col][row])
    }

    /// Extract submatrix
    #[inline]
    pub fn submat<const SR: usize, const SC: usize>(
        &self,
        row_range: Range<usize>,
        col_range: Range<usize>,
    ) -> Matrix<SR, SC> {
        debug_assert_eq!(row_range.len(), SR);
        debug_assert_eq!(col_range.len(), SC);
        debug_assert!(row_range.end <= R);
        debug_assert!(col_range.end <= C);

        let cols = array::from_fn(|col| {
            array::from_fn(|row| self.cols[col_range.start + col][row_range.start + row])
        });

        Matrix { cols }
    }

    /// Get reference to columns array
    #[inline]
    pub fn as_cols(&self) -> &[[Float; R]; C] {
        &self.cols
    }

    /// Get mutable reference to columns array
    #[inline]
    pub fn as_cols_mut(&mut self) -> &mut [[Float; R]; C] {
        &mut self.cols
    }

    /// Get columns array
    #[inline]
    pub fn to_cols(&self) -> [[Float; R]; C] {
        self.cols.clone()
    }

    /// Get raw data as flat array in row-major order
    #[inline]
    pub fn to_rows(&self) -> [[Float; C]; R] {
        let mut result = [[0.0; C]; R];

        for col in 0..C {
            for row in 0..R {
                result[row][col] = self.cols[col][row];
            }
        }

        result
    }

    /// Transpose matrix (creates new matrix with swapped dimensions)
    #[inline]
    pub fn transpose(&self) -> Matrix<C, R> {
        let cols = array::from_fn(|col| array::from_fn(|row| self.cols[row][col]));

        Matrix { cols }
    }

    /// Map function over all elements
    #[inline]
    pub fn map<F>(&self, f: F) -> Self
    where
        F: Fn(Float) -> Float,
    {
        Self {
            cols: array::from_fn(|col| array::from_fn(|row| f(self.cols[col][row]))),
        }
    }

    /// Component-wise multiplication (Hadamard product)
    #[inline]
    pub fn component_mul(&self, rhs: &Self) -> Self {
        Self {
            cols: array::from_fn(|col| {
                array::from_fn(|row| self.cols[col][row] * rhs.cols[col][row])
            }),
        }
    }

    /// Component-wise division
    #[inline]
    pub fn component_div(&self, rhs: &Self) -> Self {
        Self {
            cols: array::from_fn(|col| {
                array::from_fn(|row| self.cols[col][row] / rhs.cols[col][row])
            }),
        }
    }

    /// Check if matrix contains any NaN
    pub fn has_nan(&self) -> bool {
        self.cols.iter().any(|col| col.iter().any(|x| x.is_nan()))
    }

    /// Check if all elements are finite
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.cols
            .iter()
            .all(|col| col.iter().all(|x| x.is_finite()))
    }
}

impl<const LR: usize, const LC: usize> Matrix<LR, LC> {
    /// naive matrix multiplication
    /// compiler will produce simd code for small matrices <= ~44
    /// omeyaga slow for big matrices
    /// for `Mat2`, `Mat3`, `Mat4` use `*` operator instead
    pub fn matmul<const RC: usize>(&self, rhs: &Matrix<LC, RC>) -> Matrix<LR, RC> {
        let cols = array::from_fn(|col_idx| {
            array::from_fn(|row_idx| {
                (0..LC)
                    .map(|k| self[(row_idx, k)] * rhs[(k, col_idx)])
                    .sum()
            })
        });

        Matrix::<LR, RC>::from_cols(cols)
    }
}

// ============================================================================
// Square matrix methods
// ============================================================================

impl<const N: usize> Matrix<N, N> {
    /// Create identity matrix
    pub fn identity() -> Self {
        let mut result = Self::zero();
        for i in 0..N {
            result[(i, i)] = 1.0;
        }
        result
    }

    /// Create diagonal matrix from array
    pub fn from_diagonal(diag: [Float; N]) -> Self {
        let mut result = Self::zero();
        for i in 0..N {
            result[(i, i)] = diag[i];
        }
        result
    }

    /// Get trace (sum of diagonal elements)
    pub fn trace(&self) -> Float {
        (0..N).map(|i| self[(i, i)]).sum()
    }

    /// Check if matrix is diagonal
    pub fn is_diagonal(&self) -> bool {
        for row in 0..N {
            for col in 0..N {
                if row != col && !self[(row, col)].is_zero_eps() {
                    return false;
                }
            }
        }

        true
    }

    /// Check if matrix is symmetric
    pub fn is_symmetric(&self) -> bool {
        for row in 0..N {
            for col in (row + 1)..N {
                if !self[(row, col)].eq_eps(&self[(col, row)]) {
                    return false;
                }
            }
        }
        true
    }

    /// Check if matrix is identity
    pub fn is_identity(&self) -> bool {
        for row in 0..N {
            for col in 0..N {
                let expected = if row == col { 1.0 } else { 0.0 };

                if !self[(row, col)].eq_eps(&expected) {
                    return false;
                }
            }
        }
        true
    }

    /// Generic determinant for `Matrix<R, C>` (slower then specialized versions for `Mat2`, `Mat3`, `Mat4`)
    pub fn det_n(&self) -> Float {
        todo!("SIMD accelerated determinant for NxN matrices")
    }

    /// Generic `Matrix<R, C>` inverse (slower then specialized versions for `Mat2`, `Mat3`, `Mat4`)
    pub fn inverse_n(&self) -> Option<Self> {
        todo!("SIMD accelerated inverse for NxN matrices")
    }
}

mod impl_operations {
    use super::*;
    use opimps::{impl_ops, impl_ops_assign, impl_uni_op};
    use std::ops::*;

    impl<const R: usize, const C: usize> Index<(usize, usize)> for Matrix<R, C> {
        type Output = Float;

        #[inline]
        fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
            &self.cols[col][row]
        }
    }

    impl<const R: usize, const C: usize> IndexMut<(usize, usize)> for Matrix<R, C> {
        #[inline]
        fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
            &mut self.cols[col][row]
        }
    }

    #[impl_ops(Add)]
    fn add<const R: usize, const C: usize>(self: Matrix<R, C>, rhs: Self) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| self.cols[col][row] + rhs.cols[col][row])
        }))
    }

    #[impl_ops_assign(AddAssign)]
    fn add_assign<const R: usize, const C: usize>(self: Matrix<R, C>, rhs: Self) {
        for col in 0..C {
            for row in 0..R {
                self.cols[col][row] += rhs.cols[col][row];
            }
        }
    }

    #[impl_ops(Sub)]
    fn sub<const R: usize, const C: usize>(self: Matrix<R, C>, rhs: Self) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| self.cols[col][row] - rhs.cols[col][row])
        }))
    }

    #[impl_ops_assign(SubAssign)]
    fn sub_assign<const R: usize, const C: usize>(self: Matrix<R, C>, rhs: Self) {
        for col in 0..C {
            for row in 0..R {
                self.cols[col][row] -= rhs.cols[col][row];
            }
        }
    }

    #[impl_ops(Mul)]
    fn mul<const R: usize, const C: usize>(self: Matrix<R, C>, scalar: Float) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| self.cols[col][row] * scalar)
        }))
    }

    #[impl_ops(Mul)]
    fn mul<const R: usize, const C: usize>(self: Float, matrix: Matrix<R, C>) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| matrix.cols[col][row] * self)
        }))
    }

    #[impl_ops_assign(MulAssign)]
    fn mul_assign<const R: usize, const C: usize>(self: Matrix<R, C>, scalar: Float) {
        for col in 0..C {
            for row in 0..R {
                self.cols[col][row] *= scalar;
            }
        }
    }

    #[impl_ops(Div)]
    fn div<const R: usize, const C: usize>(self: Matrix<R, C>, scalar: Float) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| self.cols[col][row] / scalar)
        }))
    }

    #[impl_ops(Div)]
    fn div<const R: usize, const C: usize>(self: Float, matrix: Matrix<R, C>) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| matrix.cols[col][row] / self)
        }))
    }

    #[impl_ops_assign(DivAssign)]
    fn div_assign<const R: usize, const C: usize>(self: Matrix<R, C>, scalar: Float) {
        for col in 0..C {
            for row in 0..R {
                self.cols[col][row] /= scalar;
            }
        }
    }

    #[impl_uni_op(Neg)]
    fn neg<const R: usize, const C: usize>(self: Matrix<R, C>) -> Matrix<R, C> {
        Matrix::<R, C>::from_cols(array::from_fn(|col| {
            array::from_fn(|row| -self.cols[col][row])
        }))
    }
}

mod impl_default {
    use super::*;

    impl<const R: usize, const C: usize> Default for Matrix<R, C> {
        fn default() -> Self {
            Self::zero()
        }
    }
}

mod impl_eq_eps {

    use crate::core::traits::EqEps;

    use super::*;

    impl<const R: usize, const C: usize> EqEps for Matrix<R, C> {
        fn eq_eps(&self, rhs: &Self) -> bool {
            for col in 0..C {
                for row in 0..R {
                    if !self.cols[col][row].eq_eps(&rhs.cols[col][row]) {
                        return false;
                    }
                }
            }
            true
        }

        fn is_zero_eps(&self) -> bool {
            self.cols
                .iter()
                .all(|col| col.iter().all(|x| x.is_zero_eps()))
        }
    }
}

mod impl_rand {
    use super::*;

    use rand::Rng;
    use rand::distr::{Distribution, StandardUniform};

    impl<const ROW: usize, const COL: usize> Distribution<Matrix<ROW, COL>> for StandardUniform {
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Matrix<ROW, COL> {
            Matrix::<ROW, COL>::from_cols(rng.random())
        }
    }
}
