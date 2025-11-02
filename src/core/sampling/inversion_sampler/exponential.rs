use crate::core::{Float, Probpart};

impl super::InversionSampler {
    /// Computes probability density function of exponential distribution at point `x`.
    /// Used for modeling waiting times between events occurring at constant rate.
    ///
    /// Parameters:
    /// - `x`: point at which to evaluate the PDF (`x ≥ 0`)
    /// - `l`: rate parameter `λ > 0` (events per unit time)
    ///
    /// Returns probability density `λe^(-λx)`
    pub fn exponential_pdf(x: Float, l: Float) -> Float {
        l * Float::exp(-l * x)
    }

    // Step 1: Start with Exponential Distribution CDF
    // The exponential distribution has PDF: f(x) = λe^(-λx) for x ≥ 0
    // Its CDF is: F(x) = ∫₀ˣ λe^(-λt) dt = [-e^(-λt)]₀ˣ = 1 - e^(-λx)
    //
    // Step 2: Apply Inverse Transform Sampling
    // Set uniform random variable u equal to the CDF:
    // u = F(x) = 1 - e^(-λx)
    //
    // Step 3: Solve for x
    // u = 1 - e^(-λx)
    // e^(-λx) = 1 - u
    // -λx = ln(1 - u)
    // x = -ln(1 - u) / λ
    //
    // Final sampling formula: x = -ln(1 - u) / λ

    /// Given uniform random variable u in `[0, 1)` and rate parameter `l`,
    /// samples from exponential distribution.
    /// Over many calls with independent uniform `u`,
    /// returned indices are exponentially distributed.
    ///
    /// Parameters:
    /// `u`: uniform random variable in `[0, 1)`
    /// `l`: rate parameter (lambda) of exponential distribution
    ///
    /// Returns exponentially distributed sample with rate `l`.
    pub fn sample_exponential(u: Probpart, l: Float) -> Float {
        -Float::ln(1.0 - u) / l
    }

    pub fn invert_sample_exponential(x: Float, l: Float) -> Probpart {
        Probpart::new(1.0 - Float::exp(-l * x))
    }
}
