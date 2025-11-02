use crate::core::{Float, Probpart, lerp};

/// Piecewise-constant struct for 1D sampling.
///
/// Represents a tabulated, non-negative function `f(x)` over a finite interval `[min, max]`,
/// supporting efficient inverse-CDF sampling and probability evaluation.
///
/// # Invariants
/// - CDF is always valid and monotonic
/// - If total integral is zero, sampling degrades to a `uniform distribution`.
/// - Assumes the domain is continuous in `[min, max)`, with N = f.len() bins.
/// - For degenerate bins (zero weight), all sampling and PDF queries are safe.
#[derive(Debug, Clone)]
pub struct PiecewiseConst1D {
    /// Unnormalized function/PDF absolute values per bin
    pub f: Vec<Float>,

    /// Normalized cumulative distribution function `(CDF[0] = 0, CDF[f.len() + 1] = 1)`
    pub cdf: Vec<Float>,

    /// Function integral over `[min, max]`
    pub f_integ: Float,

    /// Function domain min
    pub min: Float,

    /// Function domain max
    pub max: Float,
}

impl PiecewiseConst1D {
    // Step 1: Define Piecewise Constant Function
    // Given n function samples f[0], f[1], ..., f[n-1] over interval [min, max]
    // Divide interval into n equal subintervals of width δ = (max - min)/n
    // The piecewise constant function is:
    // F(x) = f[i] for x ∈ [min + i×δ, min + (i+1)×δ), i = 0,1,...,n-1
    //
    // Step 2: Compute Cumulative Integral
    // The cumulative integral from min to point x_i = min + i×δ is:
    // I[i] = ∫[min to x_i] F(t) dt = Σ[j=0 to i-1] f[j] × δ
    //      = Σ[j=0 to i-1] f[j] × (max - min)/n
    //
    // Step 3: Normalize to Create CDF
    // Total integral: I_total = ∫[min to max] F(t) dt = Σ[j=0 to n-1] f[j] × δ
    // Normalized CDF: CDF[i] = I[i] / I_total for i = 0,1,...,n
    // where CDF[0] = 0 and CDF[n] = 1
    //
    // Step 4: Handle Degenerate Case
    // If I_total = 0 (all function values zero):
    // CDF[i] = i/n for i = 1,2,...,n
    // This creates uniform distribution over the domain
    //
    /// Constructs a piecewise constant probability distribution from function samples.
    /// Takes absolute values and normalizes to create a valid CDF for sampling.
    ///
    /// - `f`: function values sampled over the interval
    /// - `min`, `max`: domain bounds where min < max
    ///
    /// Returns `PiecewiseConst1D` with normalized `cdf`, total integral `f_integ`, and domain bounds for inverse transform sampling.
    pub fn from_func(f: &[Float], min: Float, max: Float) -> Self {
        assert!(min < max);

        let n = f.len();
        assert!(n != 0);

        let f: Vec<Float> = f.iter().map(|v| v.abs()).collect();
        let mut f_cum_integ = vec![0.0; n + 1];

        let delta = (max - min) / n as Float;

        f_cum_integ[0] = 0.0;
        for i in 1..=n {
            assert!(f[i - 1].is_finite(), "function values must be finite");
            f_cum_integ[i] = f_cum_integ[i - 1] + f[i - 1] * delta;
        }

        let f_integ = f_cum_integ[n];
        let mut cdf = f_cum_integ;

        if f_integ == 0.0 {
            // all function values zero => uniform distribution over the domain
            for i in 1..=n {
                cdf[i] = i as Float / n as Float
            }
        } else {
            for i in 1..=n {
                cdf[i] /= f_integ
            }
        }

        Self {
            f,
            cdf,
            f_integ,
            min,
            max,
        }
    }

    // Step 1: Inverse Transform Sampling Setup
    // Given uniform random variable u ∈ [0, 1) and piecewise constant function F(x)
    // over domain [min, max] divided into n equal bins of width δ = (max - min)/n
    // Each bin i has constant value f[i] and cumulative probabilities CDF[i]
    //
    // Step 2: Locate Target Bin
    // Find bin index i such that CDF[i] ≤ u < CDF[i+1]
    // This determines which constant segment of the piecewise function contains our sample
    //
    // Step 3: Remap to Local Bin Coordinate
    // Within bin i, u spans the interval [CDF[i], CDF[i+1])
    // Remap to local uniform coordinate: u_local = (u - CDF[i])/(CDF[i+1] - CDF[i])
    // When CDF[i+1] - CDF[i] = 0 (zero probability bin), u_local = 0
    //
    // Step 4: Convert to Domain Coordinate
    // Bin i spans x-interval [min + i×δ, min + (i+1)×δ)
    // Sample position within bin: x = min + (i + u_local) × δ
    //                               = min + (i + u_local) × (max - min)/n
    // Using linear interpolation: x = lerp((i + u_local)/n, min, max)
    //
    // Step 5: Compute Probability Density Function
    // PDF at x equals normalized function value: pdf = f[i]/f_integ
    // For degenerate case (f_integ = 0): pdf = 1/(max - min) (uniform density)

    /// Samples a value from the piecewise constant distribution using inverse transform sampling.
    /// Maps uniform random variable onto the normalized function distribution.
    ///
    /// - `u`: uniform random variable in `[0, 1)`
    ///
    /// Returns (`x`, `offset`, `pdf`) where `x` is sampled domain value, `offset` is bin index, `pdf` is probability density at x.
    pub fn sample(&self, u: Probpart) -> (Float, usize, Float) {
        let idx = self
            .cdf
            .partition_point(|cdf_i| cdf_i <= &u)
            .saturating_sub(1)
            .min(self.cdf.len() - 2);

        let cdf_delta = self.cdf[idx + 1] - self.cdf[idx];

        let mut u_remapped = u - self.cdf[idx];
        if cdf_delta > 0.0 {
            u_remapped /= cdf_delta;
        }

        assert!(u_remapped.is_finite());

        let pdf = if self.f_integ > 0.0 {
            self.f[idx] / self.f_integ
        } else {
            1.0 / (self.max - self.min)
        };

        let x = lerp(
            Probpart::new((idx as Float + u_remapped) / self.f.len() as Float),
            self.min,
            self.max,
        );

        (x, idx, pdf)
    }

    // As always, this is a matter of evaluating the CDF at the given position
    //
    //// Step 1: Domain to Bin Coordinate Transformation
    // Given x ∈ [min, max], normalize to [0, 1]:
    // x_norm = (x - min) / (max - min)
    //
    // Map to bin coordinates [0, n] where n = f.len():
    // c = x_norm × n = (x - min) / (max - min) × n
    //
    // Step 2: Bin Index and Fractional Position
    // idx = floor(c) ∈ [0, n-1]
    // bin_fraction = c - idx ∈ [0, 1)
    //
    // Step 3: CDF Inversion via Linear Interpolation
    // For piecewise constant functions, CDF is piecewise linear within bins.
    // The cumulative probability at position c is:
    // u = lerp(bin_fraction, cdf[idx], cdf[idx+1])
    // u = (1 - bin_fraction) × cdf[idx] + bin_fraction × cdf[idx+1]
    //
    // This transforms domain value x back to the uniform random variable u
    // that would generate x through forward sampling.

    /// Given domain value `x`, returns the cumulative probability that would sample to `x`.
    /// This inverts the sampling process - useful for warping and importance sampling.
    ///
    /// Returns `Some(u)` where `u` is the uniform random variable in `[0,1)` that maps to `x`,
    /// or `None` if `x` is outside the function domain `[min, max]`.
    ///
    pub fn invert(&self, x: Float) -> Option<Probpart> {
        if !(self.min..=self.max).contains(&x) {
            return None;
        }

        let x_norm = (x - self.min) / (self.max - self.min);
        let c = x_norm * self.f.len() as Float;

        let idx = (c.floor() as usize).clamp(0, self.f.len() - 1);
        let bin_fraction = Probpart::new(c - idx as Float);

        let u = lerp(bin_fraction, self.cdf[idx], self.cdf[idx + 1]);

        Some(Probpart::new(u))
    }
}
