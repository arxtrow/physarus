use crate::core::{Float, Point2, Probpart, lerp};

impl super::InversionSampler {
    // Step 1: Define Bilinear Interpolation on Unit Square
    // For a unit square [0,1] × [0,1] with corner weights:
    // w[0] at (0,0), w[1] at (1,0), w[2] at (0,1), w[3] at (1,1)
    //
    // The bilinear interpolation is:
    // f(x,y) = (1-x)(1-y)w[0] + x(1-y)w[1] + (1-x)y w[2] + x y w[3]
    //
    // Step 2: Compute Normalization Factor
    // For a valid PDF: ∫₀¹ ∫₀¹ pdf(x,y) dx dy = 1
    //
    // ∫₀¹ ∫₀¹ f(x,y) dx dy = ∫₀¹ ∫₀¹ [(1-x)(1-y)w[0] + x(1-y)w[1] + (1-x)y w[2] + xy w[3]] dx dy
    //
    // Each term integrates to:
    // ∫₀¹ ∫₀¹ (1-x)(1-y)w[0] dx dy = w[0] × ∫₀¹ (1-x) dx × ∫₀¹ (1-y) dy = w[0] × (1/2) × (1/2) = w[0]/4
    // ∫₀¹ ∫₀¹ x(1-y)w[1] dx dy = w[1]/4
    // ∫₀¹ ∫₀¹ (1-x)y w[2] dx dy = w[2]/4
    // ∫₀¹ ∫₀¹ xy w[3] dx dy = w[3]/4
    //
    // Total integral: (w[0] + w[1] + w[2] + w[3])/4 = w_sum/4
    //
    // Step 3: Normalize to Create PDF
    // pdf(x,y) = f(x,y) / (w_sum/4) = 4 × f(x,y) / w_sum
    //
    /// Computes bilinear probability density function at 2D point using corner weights.
    /// Returns 0.0 outside unit square, 1.0 for zero weights, normalized PDF otherwise.
    ///
    /// # Parameters
    /// * `p` - 2D point coordinates `(x, y)`
    /// * `w` - Corner weights `[bottom-left, bottom-right, top-left, top-right]`
    ///
    /// # Returns
    /// PDF value at point p
    pub fn bilinear_pdf(p: &Point2, w: &[Float; 4]) -> Float {
        let (x, y) = (p.x(), p.y());

        if !(0.0..1.0).contains(&x) || !(0.0..1.0).contains(&y) {
            return 0.0;
        }

        let w_sum: Float = w.iter().sum();

        if w_sum == 0.0 {
            return 1.0;
        }

        let fx = (1.0 - x) * (1.0 - y) * w[0]
            + x * (1.0 - y) * w[1]
            + (1.0 - x) * y * w[2]
            + x * y * w[3];

        let normalization = w_sum / 4.0;

        fx / normalization
    }

    // Step 1: Derive Marginal Distribution for y = u[1]
    // For bilinear f(x,y) = (1-x)(1-y)w[0] + x(1-y)w[1] + (1-x)y w[2] + xy w[3]
    // Marginal f_Y(y) = ∫₀¹ f(x,y) dx
    //
    // f_Y(y) = ∫₀¹ [(1-x)(1-y)w[0] + x(1-y)w[1] + (1-x)y w[2] + xy w[3]] dx
    //        = (1-y)w[0]∫₀¹(1-x)dx + (1-y)w[1]∫₀¹x dx + y w[2]∫₀¹(1-x)dx + y w[3]∫₀¹x dx
    //        = (1-y)w[0]×(1/2) + (1-y)w[1]×(1/2) + y w[2]×(1/2) + y w[3]×(1/2)
    //        = (1-y)(w[0] + w[1])/2 + y(w[2] + w[3])/2
    //        = (1/2) × [(1-y)(w[0] + w[1]) + y(w[2] + w[3])]
    //
    // Proportional form: linear distribution between (w[0] + w[1]) and (w[2] + w[3])
    // no need for (1/2), `sample_linear` handles normalization internally
    //
    // Step 2: Derive Conditional Distribution for x = u[0] given y = u[1]
    // f(x|y) = f(x,y) / f_Y(y)
    // Numerator = (1-x)[(1-y)w[0] + y w[2]] + x[(1-y)w[1] + y w[3]]
    //           = (1-x)[lerp(y, w[0], w[2])] + x[lerp(y, w[1], w[3])]
    // Proportional form: linear distribution between lerp(y, w[0], w[2]) and lerp(y, w[1], w[3])
    // no need for (/ f_Y(y)), `sample_linear` handles normalization internally
    //
    // Step 3: Two-Stage Sampling Algorithm
    // 1. Sample y from linear(w[0] + w[1], w[2] + w[3])
    // 2. Sample x from linear(lerp(y, w[0], w[2]), lerp(y, w[1], w[3]))
    ///
    /// Samples 2D point from bilinear distribution using marginal-conditional approach.
    /// First samples y coordinate, then samples x conditional on y value.
    ///
    /// # Parameters
    /// * `u` - Two independent uniform random variables [u_x, u_y]
    /// * `w` - Corner weights [bottom-left, bottom-right, top-left, top-right]
    ///
    /// # Returns
    /// Point2 sampled from bilinear distribution
    pub fn sample_bilinear(u: &[Probpart; 2], w: &[Float; 4]) -> Point2 {
        let py = Self::sample_linear(u[1], w[0] + w[1], w[2] + w[3]);
        let px = Self::sample_linear(
            u[0],
            lerp(py.into(), w[0], w[2]),
            lerp(py.into(), w[1], w[3]),
        );

        Point2::new(px, py)
    }

    /// Inverts bilinear importance sampling: given point `p` in `[0,1]²` and corner weights `w`,
    /// computes uniform variables `(u,v)` that would sample to `p` under the bilinear distribution defined by `w`.
    ///
    /// - `p`: `&Point2` - sampled point in unit square
    /// - `w`: `&[Float; 4]` - non-negative weights at corners `(0,0), (1,0), (0,1), (1,1)`
    ///
    /// Returns `(u, v)` where `u` inverts conditional x-sampling given y, `v` inverts marginal y-sampling.
    pub fn invert_sample_bilinear(p: &Point2, w: &[Float; 4]) -> (Probpart, Probpart) {
        (
            Self::invert_sample_linear(
                p.x(),
                lerp(p.y().into(), w[0], w[2]),
                lerp(p.y().into(), w[1], w[3]),
            ),
            Self::invert_sample_linear(p.y(), w[0] + w[1], w[2] + w[3]),
        )
    }
}
