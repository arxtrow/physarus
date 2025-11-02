use std::cell::Cell;

use crate::prelude::*;

#[derive(Debug, Default, Clone)]
pub struct Ray3<Space: SpaceSeal> {
    pub o: Point3<Space>,
    pub d: Unit3<Space>,
    pub time: Float,
    pub t_max: Cell<Float>,
    pub differential: Option<RayDifferential<Space>>,
}

impl<Space: SpaceSeal> Ray3<Space> {
    pub fn new(o: Point3<Space>, d: Unit3<Space>, time: Float, t_max: Float) -> Self {
        Self {
            o,
            d,
            time,
            t_max: Cell::new(t_max),
            differential: None,
        }
    }

    pub fn at(&self, t: Float) -> Point3<Space> {
        self.o + self.d * t
    }

    pub fn scale_differentials(&mut self, s: Float) {
        if let Some(ref mut diff) = self.differential {
            diff.rx_origin = self.o + (diff.rx_origin - self.o) * s;
            diff.rx_direction = self.d + (diff.rx_direction - self.d) * s;

            diff.ry_origin = self.o + (diff.ry_origin - self.o) * s;
            diff.ry_direction = self.d + (diff.ry_direction - self.d) * s;
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct RayDifferential<Space: SpaceSeal> {
    pub rx_origin: Point3<Space>,
    pub rx_direction: Direction3<Space>,

    pub ry_origin: Point3<Space>,
    pub ry_direction: Direction3<Space>,
}
