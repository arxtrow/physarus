use std::{any::Any, marker::PhantomData};

use crate::prelude::*;

pub type Vec2<Sem, Scale: ScaleSeal, Space: SpaceSeal> = Vector<2, Sem, Scale, Space>;

pub type Direction2<Space: SpaceSeal = _WorldSpace> = Vec2<_Direction, _NonUnit, Space>;
pub type Unit2<Space: SpaceSeal = _WorldSpace> = Vec2<_Direction, _Unit, Space>;

pub type Normal2<Space: SpaceSeal = _WorldSpace> = Vec2<_Normal, _Unit, Space>;

impl<Space: SpaceSeal> Normal2<Space> {
    pub fn new_normalize(x: Float, y: Float) -> Self {
        let len_sqr = x * x + y * y;

        debug_assert!(len_sqr != 0.0, "Cannot normalize zero-length vector");

        let len = len_sqr.sqrt();

        // SAFETY<Unit>: normalized by len
        unsafe { Self::raw([x / len, y / len]) }
    }

    pub fn from_vec(v1: &Vector<2, impl GeoSeal, impl ScaleSeal>) -> Self {
        let len_sqr = v1.x() * v1.x() + v1.y() * v1.y();

        debug_assert!(len_sqr != 0.0, "Cannot normalize zero-length vector");

        let len = len_sqr.sqrt();

        // SAFETY<Unit>: normalized by len
        unsafe { Self::raw([-v1.y() / len, v1.x() / len]) }
    }
}

impl<Space: SpaceSeal> Vec2<_Direction, _NonUnit, Space> {
    pub const ZERO: Self = Vec2 {
        components: [0.0, 0.0],
        phantom: PhantomData,
    };

    pub const ONE: Self = Vec2 {
        components: [1.0, 1.0],
        phantom: PhantomData,
    };
}

impl<Space: SpaceSeal> Vec2<_Direction, _Unit, Space> {
    pub const I: Self = Vec2 {
        components: [1.0, 0.0],
        phantom: PhantomData,
    };

    pub const J: Self = Vec2 {
        components: [0.0, 1.0],
        phantom: PhantomData,
    };
}

impl<Sem, Space: SpaceSeal> Vec2<Sem, _NonUnit, Space> {
    pub const fn new(x: Float, y: Float) -> Self {
        Vector {
            components: [x, y],
            phantom: PhantomData,
        }
    }
}

impl<Sem, Scale, Space> Vec2<Sem, Scale, Space>
where
    Scale: ScaleSeal,
    Space: SpaceSeal,
{
    pub fn x(&self) -> Float {
        self[0]
    }

    pub fn y(&self) -> Float {
        self[1]
    }
}

impl<SemLhs, ScaleLhs, Space> Vec2<SemLhs, ScaleLhs, Space>
where
    ScaleLhs: ScaleSeal,
    SemLhs: GeoSeal,
    Space: SpaceSeal,
{
    pub fn angle(&self, rhs: &Vec2<impl GeoSeal, impl ScaleSeal, Space>) -> Rad {
        let sin = self.x() * rhs.y() - self.y() * rhs.x(); // == |self||rhs| sin(angle)
        let cos = self.dot(rhs); // == |self||rhs| cos(angle)

        Rad::new(sin.atan2(cos).abs())
    }
}
