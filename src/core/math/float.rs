#[cfg(not(feature = "f64"))]
pub type Float = f32;

#[cfg(not(feature = "f64"))]
mod consts {
    use super::Float;

    pub const PI: Float = std::f32::consts::PI;
    pub const TAU: Float = std::f32::consts::TAU;
    pub const E: Float = std::f32::consts::E;
    pub const FRAC_PI_2: Float = std::f32::consts::FRAC_PI_2;
    pub const FRAC_PI_3: Float = std::f32::consts::FRAC_PI_3;
    pub const FRAC_PI_4: Float = std::f32::consts::FRAC_PI_4;
    pub const FRAC_PI_6: Float = std::f32::consts::FRAC_PI_6;
    pub const FRAC_PI_8: Float = std::f32::consts::FRAC_PI_8;
    pub const FRAC_1_PI: Float = std::f32::consts::FRAC_1_PI;
    pub const FRAC_2_PI: Float = std::f32::consts::FRAC_2_PI;
    pub const FRAC_2_SQRT_PI: Float = std::f32::consts::FRAC_2_SQRT_PI;
    pub const SQRT_2: Float = std::f32::consts::SQRT_2;
    pub const FRAC_1_SQRT_2: Float = std::f32::consts::FRAC_1_SQRT_2;
    pub const LN_2: Float = std::f32::consts::LN_2;
    pub const LN_10: Float = std::f32::consts::LN_10;
    pub const LOG2_E: Float = std::f32::consts::LOG2_E;
    pub const LOG10_E: Float = std::f32::consts::LOG10_E;

    pub const INFINITY: Float = f32::INFINITY;
    pub const NEG_INFINITY: Float = f32::NEG_INFINITY;
    pub const NAN: Float = f32::NAN;
    pub const EPSILON: Float = f32::EPSILON;
    pub const CLOSEST_TO_ONE: Float = 1.0_f32.next_down();
    pub const MIN: Float = f32::MIN;
    pub const MAX: Float = f32::MAX;
}

#[cfg(feature = "f64")]
pub type Float = f64;

#[cfg(feature = "f64")]
mod consts {
    use super::Float;

    pub const PI: Float = std::f64::consts::PI;
    pub const PI_OVER_2: Float = std::f64::consts::PI / 2.0;
    pub const PI_OVER_3: Float = std::f64::consts::PI / 3.0;
    pub const PI_OVER_4: Float = std::f64::consts::PI / 4.0;
    pub const TAU: Float = std::f64::consts::TAU;
    pub const E: Float = std::f64::consts::E;
    pub const FRAC_PI_2: Float = std::f64::consts::FRAC_PI_2;
    pub const FRAC_PI_3: Float = std::f64::consts::FRAC_PI_3;
    pub const FRAC_PI_4: Float = std::f64::consts::FRAC_PI_4;
    pub const FRAC_PI_6: Float = std::f64::consts::FRAC_PI_6;
    pub const FRAC_PI_8: Float = std::f64::consts::FRAC_PI_8;
    pub const FRAC_1_PI: Float = std::f64::consts::FRAC_1_PI;
    pub const FRAC_2_PI: Float = std::f64::consts::FRAC_2_PI;
    pub const SQRT_2: Float = std::f64::consts::SQRT_2;
    pub const FRAC_2_SQRT_PI: Float = std::f64::consts::FRAC_2_SQRT_PI;
    pub const LN_2: Float = std::f64::consts::LN_2;
    pub const LN_10: Float = std::f64::consts::LN_10;
    pub const LOG2_E: Float = std::f64::consts::LOG2_E;
    pub const LOG10_E: Float = std::f64::consts::LOG10_E;

    pub const INFINITY: Float = f64::INFINITY;
    pub const NEG_INFINITY: Float = f64::NEG_INFINITY;
    pub const NAN: Float = f64::NAN;
    pub const EPSILON: Float = f64::EPSILON;
    pub const CLOSEST_TO_ONE: Float = 1.0_f64.next_down();
    pub const MIN: Float = f64::MIN;
    pub const MAX: Float = f64::MAX;

    pub const ZERO: Float = 0.0;
    pub const ONE: Float = 1.0;
    pub const NEG_ONE: Float = -1.0;
}

pub use consts::*;
