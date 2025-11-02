use crate::core::{FRAC_PI_2, FRAC_PI_4, PI, Point2, Probpart};

impl super::InversionSampler {
    /// Why r = u1 is wrong (the naive mapping)
    // - In polar coordinates, area grows with radius: dA = r dr dθ.
    // - Equal steps in r don’t cover equal area: small radii cover less area per step than large radii.
    // - Making r uniform therefore crowds samples near the center and thins them near the rim./
    //
    // Step 1: Compute Joint Density p(x,y)
    // For uniform distribution over the unit disk (x² + y² ≤ 1, area π), the density is p(x,y) = 1/π inside the disk, 0 otherwise.
    //
    // Step 2: Transform to Polar Coordinates to Get p(r,θ)
    // Change variables: x = r cos θ, y = r sin θ, with Jacobian |J| = r.
    // Thus, p(r,θ) = p(x,y) ⋅ |J| = (1/π) ⋅ r, for 0 ≤ r ≤ 1, 0 ≤ θ < 2π.
    //
    // Step 3: Compute Densities
    // Marginal p(r) = ∫_0^{2π} (1/π) r dθ = (2π/π) r = 2r, for 0 ≤ r ≤ 1.
    // Conditional p(θ|r) = p(r,θ) / p(r) = [(1/π) ⋅ r] / (2r) = 1/(2π).
    //
    // Step 4: Inversion of CDF to Sample
    // For r: CDF_r(s) = ∫_0^s 2t dt = s², so inverse r = √u1
    // For θ: CDF_θ(s) = ∫_0^s 1/2π dt = 1/ 2π, so inverse θ = 2π * u2
    // Then x = r cos θ, y = r sin θ for point (x,y).
    //
    /// Samples a point uniformly from the unit disk using two independent uniform random variables `u1` and `u2`.
    /// Maps `u1` to radial distance and `u2` to angular position in polar coordinates for uniform area sampling.
    ///
    /// Although this mapping samples uniformly over the disk, it distorts areas: regions on the square may be stretched or compressed on the disk,
    /// which can reduce the effectiveness of stratified patterns. For less distortion, consider a concentric mapping from the square to the disk
    ///
    /// # Parameters
    /// - `u1`: Probpart, uniform random input for radial component.
    /// - `u2`: Probpart, uniform random input for angular component.
    ///
    /// # Returns
    /// Point2 coordinates of the sampled point in the unit disk.
    pub fn sample_disk(u1: Probpart, u2: Probpart) -> Point2 {
        let r = u1.p().sqrt();
        let theta = 2.0 * PI * u2;

        Point2::new(r * theta.cos(), r * theta.sin())
    }

    // Inverts the mapping from `sample_disk` to recover uniform samples (u₁, u₂) from a disk point p = (x, y).
    // Step 1:
    //   r = sqrt(u₁) =>  u₁ = r² = x²+y²
    //   θ = 2π*u₂    => u₂ = θ / 2π
    //
    /// Recovers the `(u₁, u₂)` uniform random variables from a disk point
    /// Given p, returns radial and angular uniforms as used by the forward mapping.
    ///
    /// # Parameters
    /// - `p`: `Point2`, point on the unit disk.
    ///
    /// # Returns
    /// `(Probpart, Probpart)`: radial and angular uniforms `(u₁, u₂)` in `[0,1)`.
    pub fn invert_sample_disk(p: Point2) -> (Probpart, Probpart) {
        let r_sqr = p.x() * p.x() + p.y() * p.y();
        let mut theta = p.y().atan2(p.x());

        if theta < 0.0 {
            theta += 2.0 * PI;
        };

        let u1 = r_sqr;
        let u2 = theta / (2.0 * PI);

        (Probpart::new(u1), Probpart::new(u2))
    }

    /// Uniform concentric square-to-disk mapping (Shirley–Chas Novak).
    /// Remaps `[0,1)²` to the unit disk with reduced distortion, preserving
    /// stratification better than polar mapping.
    ///
    /// # Core idea
    /// - First shift to `[-1,1]²`: `(v1, v2) = (2u1 - 1, 2u2 - 1)`.
    /// - The **dominant coordinate** (`|v1|` vs. `|v2|`) sets the radius.
    /// - The **ratio** of the smaller to the larger (`v2/v1` or `v1/v2`)
    ///   sets the angular offset:
    ///     - Ratio near 0 → small angular spread → thin rectangular sector.
    ///     - Ratio near ±1 → large angular spread → wide wedge-shaped sector.
    /// - This keeps mapped regions’ shapes more consistent across the disk.
    ///
    /// # Output
    /// Returns a `Point2` `(x, y)` on the unit disk (`x² + y² ≤ 1`).
    pub fn sample_disk_concentric(u1: Probpart, u2: Probpart) -> Point2 {
        let v1 = 2.0 * u1 - 1.0;
        let v2 = 2.0 * u2 - 1.0;

        if v1 == 0.0 && v2 == 0.0 {
            return Point2::ZERO;
        }

        let (r, theta) = if v1.abs() > v2.abs() {
            // Horizontal dominant
            (v1, (FRAC_PI_4) * (v2 / v1))
        } else {
            // Vertical dominant
            (v2, (FRAC_PI_2 - FRAC_PI_4 * (v1 / v2)))
        };

        Point2::new(r * theta.cos(), r * theta.sin())
    }
}
