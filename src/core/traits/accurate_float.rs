use std::simd::{Simd, StdFloat};

use crate::prelude::*;

/// Numerically stable floating-point operations for the `Float` type.
/// <https://pharr.org/matt/blog/2019/11/03/difference-of-floats.html>
pub trait AccurateFloat {
    /// Computes (self.0 × self.1) - (rhs.0 × rhs.1) with high numerical accuracy.
    fn diff_of_mul(self, rhs: Self) -> Float;

    /// Computes (self.0 × self.1) + (rhs.0 × rhs.1) using FMA.
    fn sum_of_mul(self, rhs: Self) -> Float;

    fn mul_compens(self) -> CompensatedFloat;

    fn add_compens(self) -> CompensatedFloat;
}

impl AccurateFloat for (Float, Float) {
    fn diff_of_mul(self, rhs: Self) -> Float {
        let (a, b) = self;
        let (c, d) = rhs;

        let cd: Float = c * d;
        let err = c.mul_add(-d, cd); // Capture rounding error from c*d
        let dop = a.mul_add(b, -cd); // Compute a*b - cd with FMA

        dop + err // Add back the rounding error
    }

    fn sum_of_mul(self, rhs: Self) -> Float {
        let (a, b) = self;
        let (c, d) = rhs;

        a.mul_add(b, c * d)
    }

    fn mul_compens(self) -> CompensatedFloat {
        let (a, b) = self;
        let ab = a * b;

        CompensatedFloat {
            val: ab,
            err: a.mul_add(-b, ab),
        }
    }

    fn add_compens(self) -> CompensatedFloat {
        let (a, b) = self;

        let s = a + b;
        let d = s - a;

        let a_err = a - (s - d);
        let b_err = b - d;

        CompensatedFloat {
            val: s,
            err: a_err + b_err,
        }
    }
}

#[derive(PartialEq, Debug, Clone, Copy, Default)]
pub struct CompensatedFloat {
    val: Float,
    err: Float,
}

impl From<CompensatedFloat> for f32 {
    fn from(CompensatedFloat { val, err }: CompensatedFloat) -> Self {
        val as f32 + err as f32
    }
}

impl From<CompensatedFloat> for f64 {
    fn from(CompensatedFloat { val, err }: CompensatedFloat) -> Self {
        val as f64 + err as f64
    }
}

pub trait SliceInnerProduct {
    fn inner_product(&self, rhs: &Self) -> Float;
}

impl SliceInnerProduct for [Float] {
    fn inner_product(&self, rhs: &Self) -> Float {
        debug_assert!(self.len() == rhs.len());

        let mut ip = CompensatedFloat::default();

        for (&a, &b) in self.iter().zip(rhs.iter()) {
            let mul = (a, b).mul_compens();
            let sum = (ip.val, mul.val).add_compens();

            ip.val = sum.val;
            ip.err = mul.err + sum.err;
        }

        ip.into()
    }
}

/// Extension trait for simd linear combinations
pub trait SimdLinearCombination<const LANES: usize, const N: usize> {
    /// Computes a linear combination: coeffs[0]*basis[0] + coeffs[1]*basis[1] + ...
    ///
    /// # Naming Rationale
    ///
    /// This method is called `linear_combination` rather than `dot` or `inner_product` because:
    /// - It performs **vertical** SIMD operations that preserve the vector width (returns `Simd<Float, LANES>`)
    /// - A dot product would imply a **horizontal** reduction to a scalar
    /// - The typical use case involves splatted scalar coefficients (e.g., for matrix-vector multiplication),
    ///   though the implementation accepts arbitrary SIMD vectors
    ///
    /// # Common Usage Pattern
    ///
    /// Typically used with splatted scalar coefficients for operations like matrix-vector
    /// multiplication where you compute: `v.x * col0 + v.y * col1 + v.z * col2 + v.w * col3`
    ///
    /// # Example
    /// ```
    /// // Matrix-vector multiplication pattern
    /// let coeffs = [
    ///     Simd::splat(v.x),
    ///     Simd::splat(v.y),
    ///     Simd::splat(v.z),
    /// ];
    /// let result = coeffs.linear_combination([col0, col1, col2]);
    /// ```
    fn linear_combination(self, basis: Self) -> Simd<Float, LANES>
    where
        std::simd::LaneCount<LANES>: std::simd::SupportedLaneCount,
        Self: Sized;
}

impl<const LANES: usize, const N: usize> SimdLinearCombination<LANES, N> for [Simd<Float, LANES>; N]
where
    std::simd::LaneCount<LANES>: std::simd::SupportedLaneCount,
{
    #[inline(always)]
    fn linear_combination(self, basis: Self) -> Simd<Float, LANES> {
        let coeffs = self;

        match N {
            0 => panic!("Cannot compute linear combination of 0 vectors"),

            1 => coeffs[0] * basis[0],
            2 => coeffs[0].mul_add(basis[0], coeffs[1] * basis[1]),
            3 => coeffs[0].mul_add(basis[0], coeffs[1].mul_add(basis[1], coeffs[2] * basis[2])),
            4 => coeffs[0].mul_add(
                basis[0],
                coeffs[1].mul_add(basis[1], coeffs[2].mul_add(basis[2], coeffs[3] * basis[3])),
            ),

            _ => {
                // Compiler will unroll for small N
                let mut result = coeffs[N - 1] * basis[N - 1];
                for i in (0..N - 1).rev() {
                    result = coeffs[i].mul_add(basis[i], result);
                }

                result
            }
        }
    }
}
