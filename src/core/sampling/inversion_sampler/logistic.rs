use crate::core::{Float, Probpart, lerp};

impl super::InversionSampler {
    pub fn logistic_pdf(x: Float, s: Float) -> Float {
        Float::exp(-x.abs() / s) / s * (1.0 + Float::exp(-x.abs() / s).powi(2))
    }

    pub fn sample_logistic(u: Probpart, s: Float) -> Float {
        -s * Float::ln(u / (1.0 - u))
    }

    pub fn invert_sample_logistic(x: Float, s: Float) -> Probpart {
        Probpart::new(1.0 / (1.0 + Float::exp(-x / s)))
    }

    pub fn logistic_pdf_on_range(x: Float, s: Float, a: Float, b: Float) -> Float {
        if !(a..b).contains(&x) {
            return 0.0;
        }

        let p = Self::logistic_pdf(x, s);
        let p_b = Self::invert_sample_logistic(b, s);
        let p_a = Self::invert_sample_logistic(a, s);

        p / (p_b - p_a)
    }

    pub fn sample_logistic_on_range(u: Probpart, s: Float, a: Float, b: Float) -> Float {
        let p_b = Self::invert_sample_logistic(b, s);
        let p_a = Self::invert_sample_logistic(a, s);

        let u_new = lerp(u, p_a.p(), p_b.p());

        let x = Self::sample_logistic(Probpart::new(u_new), s);

        x.clamp(a, b)
    }

    pub fn invert_sample_logistic_on_range(x: Float, s: Float, a: Float, b: Float) -> Probpart {
        let p_b = Self::invert_sample_logistic(b, s).p();
        let p_a = Self::invert_sample_logistic(a, s).p();
        let p = Self::invert_sample_logistic(x, s).p();

        Probpart::new((p - p_a) / (p_b - p_a))
    }
}
