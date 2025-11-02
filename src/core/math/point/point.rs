use std::array;
use std::marker::PhantomData;

use crate::prelude::*;

pub type Point2<Space: SpaceSeal = _WorldSpace> = Point<2, Space>;
pub type Point3<Space: SpaceSeal = _WorldSpace> = Point<3, Space>;
pub type Point4<Space: SpaceSeal = _WorldSpace> = Point<4, Space>;

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Point<const N: usize, Space: SpaceSeal = _WorldSpace> {
    pub(crate) components: [Float; N],
    state: PhantomData<(Space)>,
}

impl<const N: usize, Space: SpaceSeal> Point<N, Space> {
    pub const ORIGIN: Self = Point::new_array([0.0; N]);
    pub const ZERO: Self = Point::new_array([0.0; N]);
    pub const ONE: Self = Point::new_array([1.0; N]);
    pub const MIN: Self = Point::new_array([Float::MIN; N]);
    pub const MAX: Self = Point::new_array([Float::MAX; N]);

    pub const fn new_array(coords: [Float; N]) -> Self {
        Point {
            components: coords,
            state: PhantomData,
        }
    }

    pub fn into_direction(&self) -> Direction<N, _NonUnit, Space> {
        Direction::<N, _NonUnit, Space>::from_array_in(self.components)
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Float> {
        self.components.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Float> {
        self.components.iter_mut()
    }

    pub fn lerp(&self, other: &Point<N, Space>, t: &Self) -> Self {
        let coords = array::from_fn(|i| self[i] * (1.0 - t[i]) + other[i] * t[i]);
        Point::new_array(coords)
    }

    pub fn distance(&self, other: &Point<N, Space>) -> Float {
        self.distance_sqr(other).sqrt()
    }

    pub fn distance_sqr(&self, other: &Point<N, Space>) -> Float {
        self.iter()
            .zip(other.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum()
    }

    pub fn min(&self, other: &Point<N, Space>) -> Self {
        let coords = array::from_fn(|i| self[i].min(other[i]));

        Point::new_array(coords)
    }

    pub fn max(&self, other: &Point<N, Space>) -> Self {
        let coords = array::from_fn(|i| self[i].max(other[i]));

        Point::new_array(coords)
    }

    pub fn map<F>(&self, f: F) -> Self
    where
        F: Fn(Float) -> Float,
    {
        let coords = array::from_fn(|i| f(self[i]));
        Point::new_array(coords)
    }

    pub fn as_slice(&self) -> &[Float; N] {
        &self.components
    }

    pub fn to_array(&self) -> [Float; N] {
        self.components.clone()
    }
}

impl<Space: SpaceSeal> Point<2, Space> {
    pub const fn new(x: Float, y: Float) -> Self {
        Point::new_array([x, y])
    }

    pub fn x(&self) -> Float {
        self[0]
    }

    pub fn y(&self) -> Float {
        self[1]
    }
}

impl<Space: SpaceSeal> Point<3, Space> {
    pub const fn new(x: Float, y: Float, z: Float) -> Self {
        Point::new_array([x, y, z])
    }

    pub fn x(&self) -> Float {
        self[0]
    }

    pub fn y(&self) -> Float {
        self[1]
    }

    pub fn z(&self) -> Float {
        self[2]
    }
}

impl<Space: SpaceSeal> Point<4, Space> {
    pub const fn new(x: Float, y: Float, z: Float, w: Float) -> Self {
        Point::new_array([x, y, z, w])
    }

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

    pub fn persp_div(&self) -> Option<Point3<Space>> {
        if self.w().is_zero_eps() {
            return None;
        }

        let inv_w = 1.0 / self.w();
        Some(Point3::new(
            self.x() * inv_w,
            self.y() * inv_w,
            self.z() * inv_w,
        ))
    }
}

mod impl_from {
    use super::*;
    use crate::core::Vector;

    impl<const N: usize, VecSem, VecScale, Space: SpaceSeal>
        From<Vector<N, VecSem, VecScale, Space>> for Point<N, Space>
    {
        fn from(v: Vector<N, VecSem, VecScale, Space>) -> Self {
            Self {
                components: v.components,
                state: PhantomData,
            }
        }
    }
}

mod impl_default {
    use super::*;

    impl<const N: usize, Space: SpaceSeal> Default for Point<N, Space> {
        fn default() -> Self {
            Point::new_array([0.0; N])
        }
    }
}

mod impl_rand {
    use super::*;

    use rand::Rng;
    use rand::distr::{Distribution, StandardUniform};

    impl<const N: usize, Space: SpaceSeal> Distribution<Point<N, Space>> for StandardUniform {
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Point<N, Space> {
            Point {
                components: rng.random(),
                state: PhantomData,
            }
        }
    }
}

mod impl_operations {
    use super::*;
    use opimps::{impl_ops, impl_ops_assign};
    use std::ops::*;

    #[impl_ops(Add)]
    fn add<const N: usize, Sem, Scale: ScaleSeal, Space: SpaceSeal>(
        self: Point<N, Space>,
        rhs: Vector<N, Sem, Scale, Space>,
    ) -> Point<N, Space> {
        let coords = array::from_fn(|i| self[i] + rhs[i]);
        Point::new_array(coords)
    }

    #[impl_ops(Sub)]
    fn sub<const N: usize, Space: SpaceSeal>(
        self: Point<N, Space>,
        rhs: Point<N, Space>,
    ) -> Vector<N, _Direction, _NonUnit, Space> {
        let components = array::from_fn(|i| self[i] - rhs[i]);
        Vector::from_array_in(components)
    }

    #[impl_ops(Sub)]
    fn sub<const N: usize, Sem, Scale: ScaleSeal, Space: SpaceSeal>(
        self: Point<N, Space>,
        rhs: Vector<N, Sem, Scale, Space>,
    ) -> Point<N, Space> {
        let coords = array::from_fn(|i| self[i] - rhs[i]);
        Point::new_array(coords)
    }

    #[impl_ops(Add)]
    fn add<const N: usize, Sem, Scale: ScaleSeal, Space: SpaceSeal>(
        self: Vector<N, Sem, Scale, Space>,
        rhs: Point<N, Space>,
    ) -> Point<N, Space> {
        rhs + self
    }

    #[impl_ops_assign(AddAssign)]
    fn add_assign<const N: usize, Sem, Scale: ScaleSeal, Space: SpaceSeal>(
        self: &mut Point<N, Space>,
        rhs: Vector<N, Sem, Scale, Space>,
    ) {
        for i in 0..N {
            self[i] += rhs[i];
        }
    }

    #[impl_ops_assign(SubAssign)]
    fn sub_assign<const N: usize, Sem, Scale: ScaleSeal, Space: SpaceSeal>(
        self: &mut Point<N, Space>,
        rhs: Vector<N, Sem, Scale, Space>,
    ) {
        for i in 0..N {
            self[i] -= rhs[i];
        }
    }

    #[impl_ops(Div)]
    fn div<const N: usize, Space: SpaceSeal>(self: Point<N, Space>, rhs: Float) -> Point<N, Space> {
        let components = array::from_fn(|i| self[i] / rhs);
        Point::new_array(components)
    }

    #[impl_ops_assign(DivAssign)]
    fn div_assign<const N: usize, Space: SpaceSeal>(
        self: &mut Point<N, Space>,
        rhs: Float,
    ) -> Point<N, Space> {
        for i in 0..N {
            self[i] /= rhs;
        }
    }

    #[impl_ops(Mul)]
    fn mul<const N: usize, Space: SpaceSeal>(self: Point<N, Space>, rhs: Float) -> Point<N, Space> {
        let components = array::from_fn(|i| self[i] * rhs);
        Point::new_array(components)
    }

    #[impl_ops(Mul)]
    fn mul<const N: usize, Space: SpaceSeal>(self: Float, rhs: Point<N, Space>) -> Point<N, Space> {
        let components = array::from_fn(|i| self * rhs[i]);
        Point::new_array(components)
    }

    #[impl_ops_assign(MulAssign)]
    fn mul_assign<const N: usize, Space: SpaceSeal>(
        self: &mut Point<N, Space>,
        rhs: Float,
    ) -> Point<N, Space> {
        for i in 0..N {
            self[i] *= rhs;
        }
    }
}

mod impl_index {
    use std::ops::{
        Index, IndexMut, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
    };

    use super::*;

    impl<const N: usize, Space: SpaceSeal> IndexMut<usize> for Point<N, Space> {
        fn index_mut(&mut self, index: usize) -> &mut Self::Output {
            &mut self.components[index]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<usize> for Point<N, Space> {
        type Output = Float;

        fn index(&self, index: usize) -> &Self::Output {
            &self.components[index]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<Range<usize>> for Point<N, Space> {
        type Output = [Float];

        fn index(&self, range: Range<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> IndexMut<Range<usize>> for Point<N, Space> {
        fn index_mut(&mut self, range: Range<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<RangeFrom<usize>> for Point<N, Space> {
        type Output = [Float];

        fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> IndexMut<RangeFrom<usize>> for Point<N, Space> {
        fn index_mut(&mut self, range: RangeFrom<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<RangeTo<usize>> for Point<N, Space> {
        type Output = [Float];

        fn index(&self, range: RangeTo<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> IndexMut<RangeTo<usize>> for Point<N, Space> {
        fn index_mut(&mut self, range: RangeTo<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<RangeFull> for Point<N, Space> {
        type Output = [Float];

        fn index(&self, _range: RangeFull) -> &Self::Output {
            &self.components[..]
        }
    }

    impl<const N: usize, Space: SpaceSeal> IndexMut<RangeFull> for Point<N, Space> {
        fn index_mut(&mut self, _range: RangeFull) -> &mut Self::Output {
            &mut self.components[..]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<RangeInclusive<usize>> for Point<N, Space> {
        type Output = [Float];

        fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> IndexMut<RangeInclusive<usize>> for Point<N, Space> {
        fn index_mut(&mut self, range: RangeInclusive<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> Index<RangeToInclusive<usize>> for Point<N, Space> {
        type Output = [Float];

        fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Space: SpaceSeal> IndexMut<RangeToInclusive<usize>> for Point<N, Space> {
        fn index_mut(&mut self, range: RangeToInclusive<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }
}
