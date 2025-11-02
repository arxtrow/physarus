use crate::core::{Float, Probpart};

impl super::InversionSampler {
    pub fn tent_pdf(x: Float, r: Float) -> Probpart {
        if x.abs() >= r {
            return Probpart::ZERO;
        }

        let normalization = r * r;

        Probpart::new((r - x.abs()) / normalization)
    }

    // Step 1: Define the Tent Function
    // f(x) = r - |x| for x ∈ [-r, r], and 0 elsewhere
    // This creates a triangular distribution with peak at x = 0
    //
    // Step 2: Split into Two Halves
    // Left half: x ∈ [-r, 0], f(x) = r - |x| = r - (-x) = r + x
    // Right half: x ∈ [0, r], f(x) = r - |x| = r - x
    //
    // Step 3: Calculate Area of Each Half
    // Left area: c = ∫_{-r}^0 (r + x) dx = [rx + x²/2]_{-r}^0 = 0 - (-r² + r²/2) = r²/2
    // Right area: c = ∫_0^r (r - x) dx = [rx - x²/2]_0^r = r² - r²/2 = r²/2
    // Total area = r², so each half has probability 0.5
    //
    // Step 4: Transform Left Half to Standard Linear Distribution
    // For x ∈ [-r, 0], density: pdf_l(x) = f(x) / c = 2(r + x)/r²
    //
    // We then need to sample `X` in [−r, 0] from the conditional density pdf_l
    // A convenient trick is to make a linear change of variables that maps [-r, 0] to [0, 1]
    //
    // Substitute y = (x + r)/r so y ∈ [0, 1] (when x = -r we have y = 0, and when x = 0, y = 1)
    // and x = -r + r*y, dx/dy = r
    // pdf_l(y) = pdf_l(x) * |dx/dy| = 2(r + (-r + r*y))/r² * r = 2r*y/r² * r = 2y
    // pdf_l(y) = 2y (linear distribution from 0 to 1)
    // Inverse transform: x = -r + ry
    //
    // Step 5: Transform Right Half to Standard Linear Distribution
    // For x ∈ [0, r], density: pdf_r(x) = 2(r - x)/r²
    // Substitute y = x/r, so y ∈ [0, 1]
    // and x = ry, dx/dy = r
    // pdf_r(y) = pdf_r(x) * |dx/dy| = 2(r - ry) / r² * r  = 2(1 - y) (linear distribution from 1 to 0)
    // Inverse transform: x = ry

    /// Given uniform random variable `u` in `[0, 1)` samples from tent distribution `f(x) = r - |x|` for `x ∈ [-r, r]`, and 0 elsewhere.
    /// Over many calls with independent uniform u, returned values follow
    /// triangular distribution centered at 0 with peak density at 0 and base width 2r.
    ///
    /// Returns `x` sampled from tent distribution in range `[-r, r]`.
    pub fn sample_tent(u: Probpart, r: Float) -> Float {
        match Self::sample_discrete(&[0.5, 0.5], u) {
            (0, _, u_new) => -r + r * Self::sample_linear(u_new, 0.0, 1.0),

            (_, _, u_new) => r * Self::sample_linear(u_new, 1.0, 0.0),
        }
    }

    // Step 1: Define the Tent Distribution
    // A tent distribution with parameter r has PDF:
    // f(x) = {(r - |x|) / r², for x in [-r..r], otherwise 0}
    // or
    // f(x) = {
    // (1/r)(1 + x/r) for -r ≤ x ≤ 0
    // (1/r)(1 - x/r) for 0 ≤ x ≤ r
    // 0 otherwise
    // }
    //
    // Step 2: Decompose into Two Linear Distributions
    // Left half (-r ≤ x ≤ 0): Map x to y = -x/r
    // Then f(-ry) = (1/r)(1 - y), giving normalized PDF g_left(y) = 1 - y
    // This is a linear distribution from 1 to 0 over
    //
    // Right half (0 ≤ x ≤ r): Map x to y = x/r
    // Then f(ry) = (1/r)(1 - y), giving normalized PDF g_right(y) = 1 - y
    // This is also a linear distribution from 1 to 0 over
    //
    // Step 3: Establish the Full CDF
    // Since each half has probability 0.5:
    // - Left half maps to u ∈ [0, 0.5]
    // - Right half maps to u ∈ [0.5, 1]
    //
    // Step 4: Invert the Sampling Process
    // For x ≤ 0:
    // y = -x/r, u_half = invert_sample_linear(y, 1.0, 0.0)
    // Since left half is mirrored, u = (1 - u_half) / 2
    //
    // For x > 0:
    // y = x/r, u_half = invert_sample_linear(y, 1.0, 0.0)
    // Map to right half: u = 0.5 + u_half / 2
    //
    /// Given a sample x from tent distribution with parameter r,
    /// returns the corresponding uniform random variable that generated it.
    ///
    /// Returns uniform probability that would generate x when passed to `sample_tent`.
    pub fn invert_sample_tent(x: Float, r: Float) -> Probpart {
        let u = if x <= 0.0 {
            (1.0 - Self::invert_sample_linear(-x / r, 1.0, 0.0)) / 2.0
        } else {
            0.5 + Self::invert_sample_linear(x / r, 1.0, 0.0) / 2.0
        };

        Probpart::new(u)
    }
}
