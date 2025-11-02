use crate::prelude::*;

use super::Float;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct Probpart(Float);

impl Probpart {
    pub const ZERO: Self = Self(0.0);
    pub const ONE: Self = Self(1.0);

    pub fn new(p: Float) -> Self {
        debug_assert!(
            (0.0..=1.0).contains(&p),
            "Probability/part out of [0,1]: {p}"
        );

        Probpart(p)
    }

    pub fn p(&self) -> Float {
        self.0
    }
}

impl From<Float> for Probpart {
    fn from(value: Float) -> Self {
        Probpart::new(value)
    }
}

use rand::Rng;
use rand::distr::{Distribution, StandardUniform};

impl Distribution<Probpart> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Probpart {
        Probpart::new(rng.random())
    }
}

impl Mul<Self> for Probpart {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::new(self.0 * rhs.0)
    }
}

impl Div<Self> for Probpart {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        let p = self.0 / rhs.0;

        debug_assert!(
            (0.0..=1.0).contains(&p),
            "Probability out of [0,1]: self: {self:?}, rhs: {rhs:?}, self / rhs: {p}"
        );

        Self(p)
    }
}

impl Sub<Self> for Probpart {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let p = self.0 - rhs.0;

        debug_assert!(
            (0.0..=1.0).contains(&p),
            "Probability out of [0,1]: self: {self:?}, rhs: {rhs:?}, self - rhs: {p}"
        );

        Self(p)
    }
}

impl Add<Self> for Probpart {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let p = self.0 + rhs.0;

        debug_assert!(
            (0.0..=1.0).contains(&p),
            "Probability out of [0,1]: self: {self:?}, rhs: {rhs:?}, self + rhs: {p}"
        );

        Self(p)
    }
}

impl Sub<Float> for Probpart {
    type Output = Float;

    fn sub(self, rhs: Float) -> Float {
        self.0 - rhs
    }
}

impl Sub<Probpart> for Float {
    type Output = Float;

    fn sub(self, rhs: Probpart) -> Float {
        rhs - self
    }
}

impl Mul<Float> for Probpart {
    type Output = Float;

    fn mul(self, rhs: Float) -> Float {
        self.0 * rhs
    }
}

impl Mul<Probpart> for Float {
    type Output = Float;

    fn mul(self, rhs: Probpart) -> Float {
        rhs * self
    }
}

impl Div<Float> for Probpart {
    type Output = Float;

    fn div(self, rhs: Float) -> Float {
        self.0 / rhs
    }
}

impl Div<Probpart> for Float {
    type Output = Float;

    fn div(self, rhs: Probpart) -> Float {
        rhs / self
    }
}

impl PartialEq<Float> for Probpart {
    fn eq(&self, other: &Float) -> bool {
        self.0.eq(other)
    }
}

impl PartialEq<Probpart> for Float {
    fn eq(&self, other: &Probpart) -> bool {
        other.0.eq(self)
    }
}

impl PartialOrd<Float> for Probpart {
    fn partial_cmp(&self, other: &Float) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(other)
    }
}

impl PartialOrd<Probpart> for Float {
    fn partial_cmp(&self, other: &Probpart) -> Option<std::cmp::Ordering> {
        other.0.partial_cmp(self)
    }
}

impl EqEps for Probpart {
    fn eq_eps(&self, b: &Self) -> bool {
        self.0.eq_eps(&b.0)
    }

    fn is_zero_eps(&self) -> bool {
        self.0.is_zero_eps()
    }
}
