use std::array;

use crate::prelude::*;

pub type Bounds2<Space: SpaceSeal> = Bounds<2, Space>;
pub type Bounds3<Space: SpaceSeal> = Bounds<3, Space>;

#[derive(Debug, Copy, Clone)]
pub struct Bounds<const N: usize, Space: SpaceSeal> {
    pub p_min: Point<N, Space>,
    pub p_max: Point<N, Space>,
}

impl<const N: usize, Space: SpaceSeal> Bounds<N, Space> {
    pub fn new(p1: Point<N, Space>, p2: Point<N, Space>) -> Self {
        let p_min = p1.min(&p2);
        let p_max = p1.max(&p2);

        Self { p_min, p_max }
    }

    pub fn union(&self, rhs: &Self) -> Self {
        Self {
            p_min: self.p_min.min(&rhs.p_min),
            p_max: self.p_max.max(&rhs.p_max),
        }
    }

    pub fn union_point(&self, p: &Point<N, Space>) -> Self {
        Self {
            p_min: self.p_min.min(p),
            p_max: self.p_max.max(p),
        }
    }

    pub fn intersect(&self, rhs: &Self) -> Self {
        Self {
            p_min: self.p_min.max(&rhs.p_min),
            p_max: self.p_max.min(&rhs.p_max),
        }
    }

    pub fn overlaps(&self, rhs: &Self) -> bool {
        let self_axis_min_max = self.p_min.iter().zip(self.p_max.iter());
        let rhs_axis_min_max = rhs.p_min.iter().zip(rhs.p_max.iter());

        self_axis_min_max
            .zip(rhs_axis_min_max)
            .all(|((self_min, self_max), (rhs_min, rhs_max))| {
                self_max >= rhs_min && self_min <= rhs_max
            })
    }

    pub fn inside(&self, p: &Point<N, Space>) -> bool {
        let self_axis_min_max = self.p_min.iter().zip(self.p_max.iter());
        let p_coords = p.iter();

        self_axis_min_max
            .zip(p_coords)
            .all(|((self_min, self_max), p_coord)| (self_min..=self_max).contains(&p_coord))
    }

    pub fn inside_exclusive(&self, p: &Point<N, Space>) -> bool {
        let self_axis_min_max = self.p_min.iter().zip(self.p_max.iter());
        let p_coords = p.iter();

        self_axis_min_max
            .zip(p_coords)
            .all(|((self_min, self_max), p_coord)| (self_min..self_max).contains(&p_coord))
    }

    pub fn distance_sqr(&self, p: &Point<N, Space>) -> Float {
        let self_axis_min_max = self.p_min.iter().zip(self.p_max.iter());
        let p_coords = p.iter();

        self_axis_min_max
            .zip(p_coords)
            .map(|((self_min, self_max), p_coord)| {
                (self_min - p_coord).max(p_coord - self_max).max(0.0)
            })
            .map(|d| d * d)
            .sum()
    }

    pub fn distance(&self, p: &Point<N, Space>) -> Float {
        self.distance_sqr(p).sqrt()
    }

    pub fn expand(&self, delta: Float) -> Self {
        Self {
            p_min: self.p_min - Direction::from_array_in([delta; N]),
            p_max: self.p_max + Direction::from_array_in([delta; N]),
        }
    }

    pub fn diagonal(&self) -> Direction<N, _NonUnit, Space> {
        self.p_max - self.p_min
    }

    pub fn max_dimension(&self) -> usize {
        let d = self.diagonal();

        (0..N)
            .max_by(|&i, &j| d[i].partial_cmp(&d[j]).unwrap())
            .unwrap_or(0)
    }

    pub fn lerp(&self, t: &Point<N, Space>) -> Point<N, Space> {
        self.p_min.lerp(&self.p_max, t)
    }

    pub fn offset(&self, p: &Point<N, Space>) -> Direction<N, _NonUnit, Space> {
        let o = p - self.p_min;

        let components = array::from_fn(|i| {
            if self.p_max[i] > self.p_min[i] {
                o[i] / (self.p_max[i] - self.p_min[i])
            } else {
                o[i]
            }
        });

        Direction::from_array_in(components)
    }

    /// Returns: `(center, radius)`
    pub fn bounding_sphere(&self) -> (Point<N, Space>, Float) {
        let center = self.p_min + self.diagonal() / 2.0;
        let radius = self.p_max.distance(&center);

        (center, radius)
    }

    pub fn is_empty(&self) -> bool {
        self.p_min
            .iter()
            .zip(self.p_max.iter())
            .any(|(min, max)| min >= max)
    }

    pub fn is_degenerate(&self) -> bool {
        self.p_min
            .iter()
            .zip(self.p_max.iter())
            .any(|(min, max)| min > max)
    }
}

impl<Space: SpaceSeal> Bounds<3, Space> {
    pub fn surface_area(&self) -> Float {
        let d = self.diagonal();

        let bottom_area = d.x() * d.y();
        let left_area = d.x() * d.z();
        let back_area = d.y() * d.z();

        2.0 * (bottom_area + left_area + back_area)
    }

    pub fn volume(&self) -> Float {
        let d = self.diagonal();

        d.x() * d.y() * d.z()
    }

    pub fn corner(&self, idx: u8) -> Point3<Space> {
        debug_assert!(idx < 8);

        let p0 = self.p_min;
        let p1 = self.p_max;

        let x = if idx & 0b001 != 0 { p1.x() } else { p0.x() };
        let y = if idx & 0b010 != 0 { p1.y() } else { p0.y() };
        let z = if idx & 0b100 != 0 { p1.z() } else { p0.z() };

        Point3::new(x, y, z)
    }
}

impl<Space: SpaceSeal> Bounds<2, Space> {
    pub fn area(&self) -> Float {
        let d = self.diagonal();

        d.x() * d.y()
    }

    pub fn corner(&self, idx: u8) -> Point2<Space> {
        debug_assert!(idx < 4);

        let p0 = self.p_min;
        let p1 = self.p_max;

        let x = if idx & 0b01 != 0 { p1.x() } else { p0.x() };
        let y = if idx & 0b10 != 0 { p1.y() } else { p0.y() };

        Point2::new(x, y)
    }
}

impl<const N: usize, Space: SpaceSeal> Default for Bounds<N, Space> {
    fn default() -> Self {
        Self {
            p_min: Point::MAX,
            p_max: Point::MIN,
        }
    }
}

impl<const N: usize, Space: SpaceSeal> From<Point<N, Space>> for Bounds<N, Space> {
    fn from(p: Point<N, Space>) -> Self {
        Self { p_min: p, p_max: p }
    }
}

mod impl_rand {
    use super::*;

    use rand::Rng;
    use rand::distr::{Distribution, StandardUniform};

    impl<const N: usize, Space: SpaceSeal> Distribution<Bounds<N, Space>> for StandardUniform {
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Bounds<N, Space> {
            Bounds::new(rng.random(), rng.random())
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::{Bounds3, Point3};

    #[test]
    fn debug_test() {
        let bb1 = Bounds3::new(rand::random::<Point3>(), rand::random::<Point3>());
        let bb2 = Bounds3::new(rand::random::<Point3>(), rand::random::<Point3>());
        dbg!(bb1.overlaps(&bb2));
    }
}
