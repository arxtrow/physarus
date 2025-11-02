use crate::prelude::*;

#[inline]
pub fn lerp(t: Probpart, a: Float, b: Float) -> Float {
    (Probpart::ONE - t) * a + t * b
}

#[inline]
pub fn smooth_step(x: Float, a: Float, b: Float) -> Float {
    if a == b {
        return if x < a { 0.0 } else { 1.0 };
    }

    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);

    t * t * (3.0 - 2.0 * t)
}
