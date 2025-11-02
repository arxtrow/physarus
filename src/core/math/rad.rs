use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Rad(pub Float);

impl Rad {
    pub const ZERO: Self = Self(0.0);
    pub const FRAC_PI_2: Self = Self(FRAC_PI_2);
    pub const PI: Self = Self(PI);
    pub const TAU: Self = Self(TAU);

    pub fn new(r: Float) -> Self {
        Rad(r)
    }

    pub fn normalize(self) -> Self {
        Self(self.0.rem_euclid(TAU))
    }

    pub fn normalize_signed(self) -> Self {
        let Rad(a) = self.normalize();

        if a > PI { Rad(a - TAU) } else { Rad(a) }
    }

    pub fn from_degrees(deg: Float) -> Self {
        Rad(deg / 360.0 * TAU)
    }

    pub fn to_degrees(self) -> Float {
        self.0 / TAU * 360.0
    }

    pub fn eq_mod_tau(self, other: Rad) -> bool {
        let Rad(a) = self;
        let Rad(b) = other;

        let d = (a - b).rem_euclid(TAU);

        d.is_zero_eps() || d.eq_eps(&TAU)
    }

    pub fn sin_cos(self) -> (Float, Float) {
        self.0.sin_cos()
    }

    pub fn sin(self) -> Float {
        self.0.sin()
    }

    pub fn cos(self) -> Float {
        self.0.cos()
    }

    pub fn tan(self) -> Float {
        self.0.tan()
    }

    pub fn asin(self) -> Self {
        Self(self.0.asin())
    }

    pub fn acos(self) -> Self {
        Self(self.0.acos())
    }

    pub fn atan(self) -> Self {
        Self(self.0.atan())
    }

    pub fn atan2(self, rhs: Rad) -> Rad {
        Self(self.0.atan2(rhs.0))
    }

    pub fn abs(&self) -> Self {
        Self(self.0.abs())
    }
}

impl Add for Rad {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Rad {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Mul<Float> for Rad {
    type Output = Self;

    fn mul(self, rhs: Float) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<Rad> for Float {
    type Output = Rad;

    fn mul(self, rhs: Rad) -> Self::Output {
        Rad(self * rhs.0)
    }
}

impl Div<Float> for Rad {
    type Output = Self;

    fn div(self, rhs: Float) -> Self::Output {
        Self(self.0 / rhs)
    }
}

impl Neg for Rad {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}
