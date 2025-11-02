use std::marker::PhantomData;

use crate::prelude::*;

pub type Vec4<Sem, Scale: ScaleSeal, Space: SpaceSeal> = Vector<4, Sem, Scale, Space>;

pub type Direction4<Space: SpaceSeal = _WorldSpace> = Vec4<_Direction, _NonUnit, Space>;

pub type Unit4<Space: SpaceSeal = _WorldSpace> = Vec4<_Direction, _Unit, Space>;

impl Direction4 {
    pub const ZERO: Self = Vec4 {
        components: [0.0, 0.0, 0.0, 0.0],
        phantom: PhantomData,
    };

    pub const ONE: Self = Vec4 {
        components: [1.0, 1.0, 1.0, 1.0],
        phantom: PhantomData,
    };
}

impl Unit4 {
    pub const I: Self = Vec4 {
        components: [1.0, 0.0, 0.0, 0.0],
        phantom: PhantomData,
    };

    pub const J: Self = Vec4 {
        components: [0.0, 1.0, 0.0, 0.0],
        phantom: PhantomData,
    };

    pub const K: Self = Vec4 {
        components: [0.0, 0.0, 1.0, 0.0],
        phantom: PhantomData,
    };

    pub const W: Self = Vec4 {
        components: [0.0, 0.0, 0.0, 1.0],
        phantom: PhantomData,
    };
}

impl<Sem, Space: SpaceSeal> Vec4<Sem, _NonUnit, Space> {
    pub const fn new(x: Float, y: Float, z: Float, w: Float) -> Self {
        Vector {
            components: [x, y, z, w],
            phantom: PhantomData,
        }
    }
}

impl<Sem, Scale, Space: SpaceSeal> Vec4<Sem, Scale, Space>
where
    Scale: ScaleSeal,
{
    pub fn x(&self) -> Float {
        self[0]
    }

    pub fn y(&self) -> Float {
        self[1]
    }

    pub fn z(&self) -> Float {
        self[2]
    }

    pub fn w(&self) -> Float {
        self[3]
    }
}
