use crate::core::{Float, Probpart, newton_bisection, smooth_step};

impl super::InversionSampler {
    pub fn smooth_step_pdf(x: Float, a: Float, b: Float) -> Float {
        if !(a..b).contains(&x) {
            return 0.0;
        }

        let normalization = 2.0 / (b - a);

        let fx = smooth_step(x, a, b);

        normalization * fx
    }

    // Step 1: Define the Smooth Step Function
    // s(t) = t²(3 - 2t) = 3t² - 2t³ for t ∈ [0,1]
    // This is a smooth interpolation curve with zero derivatives at endpoints
    //
    // Step 2: Normalize to Create a Probability Density Function (PDF)
    // ∫₀¹ s(t) dt = ∫₀¹ (3t² - 2t³) dt = [t³ - t⁴/2]₀¹ = 1 - 1/2 = 1/2
    //
    // To normalize: pdf(t) = s(t) / 1/2 = 2(3t² - 2t³) = 6t² - 4t³
    //
    // Step 3: Derive the Cumulative Distribution Function (CDF)
    // CDF(t) = ∫₀ᵗ pdf(u) du = ∫₀ᵗ (6u² - 4u³) du = [2u³ - u⁴]₀ᵗ = 2t³ - t⁴
    //
    // For x ∈ [a,b] with t = (x-a)/(b-a):
    // CDF(x) = 2t³ - t⁴ where t = (x-a)/(b-a)
    //
    // Step 4: Invert CDF for Sampling
    // Given u ∈ [0,1), solve: CDF(x) = u
    // 2t³ - t⁴ = u
    // This quartic equation requires numerical root-finding methods
    // Using Newton bisection find x for CDF(x) - u = 0
    //
    /// Given uniform random variable `u` in `[0, 1)`, samples from smooth step distribution.
    /// The distribution has density `pdf(x) ∝ 3t² - 2t³` where `t = (x-a)/(b-a)`.
    ///
    /// - `u`: uniform random variable in `[0, 1)`
    /// - `a`: lower bound of distribution
    /// - `b`: upper bound of distribution
    ///
    /// Returns x sampled from smooth step distribution over `[a, b]`.
    pub fn sample_smooth_step(u: Probpart, a: Float, b: Float) -> Float {
        let cdf_minus_u = |x: Float| {
            let t = (x - a) / (b - a);
            let p = 2.0 * t.powi(3) - t.powi(4);

            let dp = Self::smooth_step_pdf(x, a, b);

            (p - u, dp)
        };

        newton_bisection(cdf_minus_u, a, b)
    }

    pub fn invert_sample_smooth_step(x: Float, a: Float, b: Float) -> Probpart {
        let t = (x - a) / (b - a);
        let p = |v: Float| 2.0 * v.powi(3) - v.powi(4);

        Probpart::new((p(t) - p(a)) / (p(b) - p(a)))
    }
}

#[cfg(test)]
mod tests {
    use crate::core::{Float, InversionSampler, Probpart, lerp, smooth_step};

    use plotters::prelude::*;

    #[test]
    fn plot_smooth_step_samples() {
        let a = 3.0_f64;
        let b = 8.0_f64;
        let samples = 1000;

        let mut points: Vec<(f64, f64)> = Vec::with_capacity(samples);

        for i in 0..samples {
            let u: Probpart = rand::random();
            let x = InversionSampler::sample_smooth_step(u, a as Float, b as Float);
            points.push((
                x as f64,
                smooth_step(x as Float, a as Float, b as Float) as f64,
            ));
        }

        let root = BitMapBackend::new("smooth_step_samples.png", (800, 600)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .caption("Smooth-Step Samples", ("sans-serif", 40))
            .margin(10)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0.0..50.0 as f64, 0.0..2.0)
            .unwrap();

        chart.configure_mesh().draw().unwrap();

        chart
            .draw_series(
                points
                    .iter()
                    .map(|(idx, val)| Circle::new((*idx, *val), 2, BLUE.filled())),
            )
            .unwrap();

        println!("Plot saved to smooth_step_samples.png");
    }
}
