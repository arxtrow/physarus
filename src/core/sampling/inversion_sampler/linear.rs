use crate::core::{CLOSEST_TO_ONE, Float, Probpart, lerp};

impl super::InversionSampler {
    fn linear_pdf(t: Probpart, a: Float, b: Float) -> Float {
        let normalization = (a + b) / 2.0;
        let fx = lerp(t, a, b);

        fx / normalization
    }

    // Step 1: Normalize to Create a Probability Density Function (PDF)
    // f(x) = (1-x)a + bx = a - ax + bx = a + (b-a)x
    // ∫₀¹ f(x) dx = ∫₀¹ [a + (b-a)x] dx
    //             = [ax + (b-a)x²/2]₀¹
    //             = a + (b-a)/2
    //             = (a + b)/2
    //
    // pdf(x) = f(x) / ((a+b)/2) = 2f(x)/(a+b) = 2[a + (b-a)x]/(a+b)
    //
    // Step 2: Derive the CDF;
    // CDF(x) = ∫₀ˣ pdf(t) dt = (2/(a+b)) ∫₀ˣ [a + (b-a)t] dt
    //                        = (2/(a+b)) [at + (b-a)t²/2]₀ˣ
    //                        = (2/(a+b)) [ax + (b-a)x²/2]
    //                        = 2ax/(a+b) + (b-a)x²/(a+b)
    //                        = x[2a + (b-a)x]/(a+b)
    //
    // Step 3: Invert the CDF for Sampling
    // u = x[2a + (b-a)x]/(a+b)
    // u(a+b) = x[2a + (b-a)x] = 2ax + (b-a)x²
    // (b-a)x² + 2ax - u(b+a) = 0
    //
    // Step 4: Apply Quadratic Formula x = (-B ± √(B² - 4AC))/(2A)
    // x = [-2a ± √(4a² + 4u(b-a)(b+a))] / [2(b-a)]
    // taking x >= 0
    // x = [-a + √(a² + u(b-a)(b+a))] / (b-a)
    //   = [-a + √(a² + u(b² - a²))] / (b-a)
    // in this form, when a = b, result is a division by 0
    //
    // Step 5: Transform to Code Formula
    // x = [-a + √(a² + u(b² - a²))] × [a + √(a² + u(b² - a²))]
    //      ────────────────────────────────────────────────────────
    //                         (b-a) × [a + √(a² + u(b² - a²))]
    //
    // The numerator becomes:
    // [-a + √(a² + u(b² - a²))] × [a + √(a² + u(b² - a²))]
    //               = [√(a² + u(b² - a²))]² - a²
    //               = a² + u(b² - a²) - a²
    //               = u(b² - a²)
    //               = u(b-a)(b+a)
    //
    // x = u(b-a)(b+a) / [(b-a)[a + √(a² + u(b² - a²))]]
    // x = u(a+b) / [a + √(a² + u(b² - a²))]
    //
    // lerp(u, a², b²) = (1 - u)a² + ub² = a² + u(b² - a²)
    //
    // x = u * (a + b) / (a + sqrt(lerp(u, a * a, b * b)))

    /// Given uniform random variable `u` in `[0, 1)`
    /// samples an index from continuous linear distribution
    /// Over many calls with independent uniform `u`,
    /// returned indices are linearly distributed between `a` and `b`.
    ///
    /// Returns `x` that by mapping `u` onto linear distribution from `a` to `b`.
    pub fn sample_linear(u: Probpart, a: Float, b: Float) -> Float {
        debug_assert!(a >= 0.0 && b >= 0.0);

        match (u, a, b) {
            // avoids division by zero below
            (Probpart::ZERO, 0.0, _) => 0.0,

            (u, a, b) => {
                let x = u * (a + b) / (a + Float::sqrt(lerp(u, a * a, b * b)));

                x.min(CLOSEST_TO_ONE)
            }
        }
    }

    // CDF(x) from Step 2
    //
    /// Returns cumulative probability that x represents in a linear distribution between `a` and `b`.
    pub fn invert_sample_linear(x: Float, a: Float, b: Float) -> Probpart {
        debug_assert!(a <= b && (a..b).contains(&x));

        Probpart::new(x * (a * (2.0 - x) + b * x) / (a + b))
    }
}
