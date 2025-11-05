use std::array;
use std::marker::PhantomData;

use crate::prelude::*;

pub type Direction<const N: usize, Scale: ScaleSeal = _NonUnit, Space: SpaceSeal = _WorldSpace> =
    Vector<N, _Direction, Scale, Space>;

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Vector<const N: usize, Sem, Scale, Space = _WorldSpace> {
    pub(in crate::core) components: [Float; N],
    pub(in crate::core) phantom: PhantomData<([(); N], Sem, Scale, Space)>,
}

impl<const N: usize, SemLhs> Vector<N, SemLhs, _NonUnit, _WorldSpace> {
    pub const fn from_array(components: [Float; N]) -> Self {
        Vector {
            components,
            phantom: PhantomData,
        }
    }
}

impl<const N: usize, SemLhs, Space: SpaceSeal> Vector<N, SemLhs, _NonUnit, Space> {
    pub const fn from_array_in(components: [Float; N]) -> Self {
        Vector {
            components,
            phantom: PhantomData,
        }
    }
}

impl<const N: usize, SemLhs, ScaleLhs, Space> Vector<N, SemLhs, ScaleLhs, Space>
where
    ScaleLhs: ScaleSeal,
    Space: SpaceSeal,
{
    pub const unsafe fn raw(components: [Float; N]) -> Self {
        Vector {
            components,
            phantom: PhantomData,
        }
    }

    pub const fn in_space<To: SpaceSeal>(&self) -> Vector<N, SemLhs, ScaleLhs, To> {
        Vector {
            components: self.components,
            phantom: PhantomData,
        }
    }

    pub fn as_slice(&self) -> &[Float; N] {
        &self.components
    }

    pub fn to_array(&self) -> [Float; N] {
        self.components.clone()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Float> {
        self.components.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Float> {
        self.components.iter_mut()
    }

    pub fn abs(&self) -> Vector<N, SemLhs, ScaleLhs> {
        let res = array::from_fn(|i| self[i].abs());

        Vector {
            components: res,
            phantom: PhantomData,
        }
    }

    pub fn has_nan(&self) -> bool {
        self.components.iter().any(|c| c.is_nan())
    }

    pub fn permute(&self, indices: [usize; N]) -> Vector<N, SemLhs, ScaleLhs, Space> {
        let components = array::from_fn(|i| self[indices[i]]);

        Vector {
            components,
            phantom: PhantomData,
        }
    }

    pub fn is_finite(&self) -> bool {
        self.components.iter().all(|c| c.is_finite())
    }

    pub fn is_zero(&self) -> bool {
        self.components.iter().all(|c| c.is_zero_eps())
    }

    pub fn clamp(&self, min: Float, max: Float) -> Vector<N, SemLhs, _NonUnit, Space> {
        let components = array::from_fn(|i| self[i].clamp(min, max));

        unsafe { Vector::raw(components) }
    }

    pub fn lerp<SemRhs, ScaleRhs, SemT, ScaleT>(
        &self,
        other: &Vector<N, SemRhs, ScaleRhs, Space>,
        t: &Vector<N, SemT, ScaleT, Space>,
    ) -> Vector<N, SemLhs, _NonUnit, Space>
    where
        ScaleRhs: ScaleSeal,
        ScaleT: ScaleSeal,
    {
        let res = array::from_fn(|i| (1.0 - t[i]) * self[i] + other[i] * t[i]);

        unsafe { Vector::raw(res) }
    }
}

impl<const N: usize, Sem, Scale, Space> Vector<N, Sem, Scale, Space>
where
    Scale: ScaleSeal,
    Space: SpaceSeal,
{
    pub fn len_sqr(&self) -> Float {
        if Scale::IS_UNIT {
            return 1.0;
        }

        self.iter().map(|val| val * val).sum()
    }

    pub fn len(&self) -> Float {
        if Scale::IS_UNIT {
            return 1.0;
        }

        self.len_sqr().sqrt()
    }
}

impl<const N: usize, SemLhs, ScaleLhs, Space> Vector<N, SemLhs, ScaleLhs, Space>
where
    SemLhs: GeoSeal,
    ScaleLhs: ScaleSeal,
    Space: SpaceSeal,
{
    pub fn is_parallel(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> bool {
        let len_self = self.len();
        let len_rhs = rhs.len();

        if len_self.is_zero_eps() || len_rhs.is_zero_eps() {
            return false;
        }

        self.cos_angle(rhs).eq_eps(&1.0)
    }

    pub fn is_perpendicular(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> bool {
        self.dot(rhs).abs().is_zero_eps()
    }

    pub fn cos_angle(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> Float {
        self.dot(rhs) / (self.len() * rhs.len())
    }

    pub fn dot(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> Float {
        self.iter()
            .zip(rhs.iter())
            .map(|(lhs, rhs)| lhs * rhs)
            .sum()
    }

    pub fn distance(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> Float {
        (self - rhs).len()
    }

    pub fn distance_sqr(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> Float {
        (self - rhs).len_sqr()
    }
}

impl<const N: usize, SemLhs, Space> Vector<N, SemLhs, _NonUnit, Space>
where
    SemLhs: GeoSeal,
    Space: SpaceSeal,
{
    pub fn angle_n(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> Rad {
        // By law of cosines:
        // |self * |rhs| - rhs * |self||^2 = |self|^2 * |rhs|^2| + |rhs|^2*|self|^2 - 2|self|^2|rhs|^2*cos(a)
        // |self * |rhs| + rhs * |self||^2 = |self|^2 * |rhs|^2| + |rhs|^2*|self|^2 + 2|self|^2|rhs|^2*cos(a)
        //
        // |self * |rhs| - rhs * |self||^2 = (1 - cos(a)) * (2*self|^2 * |rhs|^2|)
        // |self * |rhs| + rhs * |self||^2 = (1 + cos(a)) * (2*self|^2 * |rhs|^2|)
        //
        // tan^2(a/2) = |1 - 2cosa| / |1 + 2cos(a)| = |self * |rhs| - rhs * |self||^2 / |self * |rhs| + rhs * |self||^2
        // => tan(a/2) = |self * |rhs| - rhs * |self|| / |self * |rhs| - rhs * |self||
        // => a = 2.0 * atan2(|self * |rhs| - rhs * |self||, |self * |rhs| - rhs * |self||)
        //
        let len_self = self.len();
        let len_rhs = rhs.len();

        Rad::new(
            2.0 * ((len_rhs * self - len_self * rhs).len())
                .atan2((len_rhs * self + len_self * rhs).len()),
        )
    }

    /// project `rhs` on `self`
    pub fn project(
        &self,
        rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Direction<N, _NonUnit, Space> {
        // project rhs on self = (rhs · (self / |self|)) × (self / |self|) = [(rhs · self) / |self|²] × self
        let len_sqr = self.len_sqr();

        debug_assert!(!len_sqr.is_zero_eps());

        (rhs.dot(self) / len_sqr) * self
    }

    /// rejection of `rhs` from `self`
    pub fn reject(
        &self,
        rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Direction<N, _NonUnit, Space> {
        rhs - self.project(rhs)
    }

    pub fn gram_schmidt(
        &self,
        rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Direction<N, _Unit, Space> {
        self.reject(rhs).normalize()
    }

    /// reflect `rhs` over `self`
    pub fn reflect<ScaleRhs: ScaleSeal>(
        &self,
        rhs: &Vector<N, impl GeoSeal, ScaleRhs, Space>,
    ) -> Direction<N, ScaleRhs, Space> {
        let reflected = rhs - 2.0 * self.reject(rhs);

        // SAFETY<Scale>: reflection preserve scale
        unsafe { Direction::raw(reflected.components) }
    }

    pub fn normalize(&self) -> Vector<N, SemLhs, _Unit, Space> {
        let normalized = self / self.len();

        // SAFETY<Unit>: normalized vector
        unsafe { Vector::raw(normalized.components) }
    }

    pub fn try_normalize(&self) -> Option<Vector<N, SemLhs, _Unit, Space>> {
        let len_sqr = self.len_sqr();

        if len_sqr == 0.0 {
            return None;
        }

        let normalized = self / len_sqr.sqrt();

        // SAFETY<Unit>: normalized vector
        unsafe { Some(Vector::raw(normalized.components)) }
    }
}

impl<const N: usize, SemLhs, Space> Vector<N, SemLhs, _Unit, Space>
where
    SemLhs: GeoSeal,
    Space: SpaceSeal,
{
    pub fn angle_n(&self, rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>) -> Rad {
        // |self| = 1 = |rhs|
        //
        // By law of cosines:
        // |self - rhs|^2 = |self|^2 + |rhs|^2 - 2|self||rhs|cos(a) = 1 - 2cos(a)
        // |self + rhs|^2 = |self|^2 + |rhs|^2 + 2|self||rhs|cos(a) = 1 + 2cos(a)
        //
        // tan^2(a/2) = |1 - 2cos(a)| / |1 + 2cos(a)| = |self - rhs|^2 / |self + rhs|^2
        // => tan(a/2) = |self - hrs| / |self + rhs|
        // => a = 2.0 * atan2(|self - rhs|, |self + rhs|)
        Rad::new(2.0 * (self - rhs).len().atan2((self + rhs).len()))
    }

    /// project `rhs` onto `self`
    pub fn project(
        &self,
        rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Direction<N, _NonUnit, Space> {
        // project rhs on self = (rhs · (self / |self|)) × (self / |self|) = [(rhs · self) / |self|²] × self
        // |self| == 1 => project rhs on self = (rhs · self) × self
        rhs.dot(self) * self
    }

    /// rejection of `rhs` from `self`
    pub fn reject(
        &self,
        rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Direction<N, _NonUnit, Space> {
        rhs - self.project(rhs)
    }

    pub fn gram_schmidt(
        &self,
        rhs: &Vector<N, impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Direction<N, _Unit, Space> {
        self.reject(rhs).normalize()
    }

    /// reflect `rhs` over `self`
    pub fn reflect<ScaleRhs: ScaleSeal>(
        &self,
        rhs: &Vector<N, impl GeoSeal, ScaleRhs, Space>,
    ) -> Direction<N, ScaleRhs, Space> {
        let reflected = rhs - 2.0 * self.reject(rhs);

        // SAFETY<Scale>: reflection preserve scale
        unsafe { Direction::raw(reflected.components) }
    }
}

mod impl_eq_eps {
    use super::*;

    impl<const N: usize, SemLhs, ScaleLhs, Space> EqEps for Vector<N, SemLhs, ScaleLhs, Space>
    where
        ScaleLhs: ScaleSeal,
        Space: SpaceSeal,
    {
        fn eq_eps(&self, rhs: &Self) -> bool {
            self.components
                .iter()
                .zip(rhs.components.iter())
                .all(|(l, r)| l.eq_eps(r))
        }

        fn is_zero_eps(&self) -> bool {
            self.components.iter().all(|c| c.is_zero_eps())
        }
    }
}

mod impl_operations {
    use super::*;
    use std::ops::*;

    use opimps::{impl_ops, impl_ops_assign, impl_ops_lprim, impl_ops_rprim, impl_uni_ops};

    #[impl_ops(Add)]
    fn add<const N: usize, SemLhs, SemRhs, ScaleLhs, ScaleRhs, Space>(
        self: Vector<N, SemLhs, ScaleLhs, Space>,
        rhs: Vector<N, SemRhs, ScaleRhs, Space>,
    ) -> Direction<N, _NonUnit, Space>
    where
        ScaleLhs: ScaleSeal,
        ScaleRhs: ScaleSeal,
        Space: SpaceSeal,
    {
        let res = array::from_fn(|i| self[i] + rhs[i]);
        Vector::from_array_in(res)
    }

    #[impl_ops(Sub)]
    fn sub<const N: usize, SemLhs, SemRhs, ScaleLhs, ScaleRhs, Space>(
        self: Vector<N, SemLhs, ScaleLhs, Space>,
        rhs: Vector<N, SemRhs, ScaleRhs, Space>,
    ) -> Direction<N, _NonUnit, Space>
    where
        ScaleLhs: ScaleSeal,
        ScaleRhs: ScaleSeal,
        Space: SpaceSeal,
    {
        let res = array::from_fn(|i| self[i] - rhs[i]);
        Vector::from_array_in(res)
    }

    #[impl_ops_rprim(Mul)]
    fn mul<const N: usize, Sem, Scale, Space>(
        self: Vector<N, Sem, Scale, Space>,
        rhs: Float,
    ) -> Direction<N, _NonUnit, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        let res = array::from_fn(|i| self[i] * rhs);
        Vector::from_array_in(res)
    }

    #[impl_ops_lprim(Mul)]
    fn mul<const N: usize, Sem, Scale, Space>(
        self: Float,
        rhs: Vector<N, Sem, Scale, Space>,
    ) -> Direction<N, _NonUnit, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        let res = array::from_fn(|i| rhs[i] * self);
        Vector::from_array_in(res)
    }

    #[impl_ops_rprim(Div)]
    fn div<const N: usize, Sem, Scale, Space>(
        self: Vector<N, Sem, Scale, Space>,
        rhs: Float,
    ) -> Direction<N, _NonUnit, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        let res = array::from_fn(|i| self[i] / rhs);
        Vector::from_array_in(res)
    }

    #[impl_uni_ops(Neg)]
    fn neg<const N: usize, Sem, Scale, Space>(
        self: Vector<N, Sem, Scale, Space>,
    ) -> Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        let components = array::from_fn(|i| -self[i]);

        Vector {
            components,
            phantom: PhantomData,
        }
    }

    #[impl_ops_assign(AddAssign)]
    fn add_assign<const N: usize, Sem, ScaleRhs, Space>(
        self: &mut Vector<N, Sem, _NonUnit, Space>,
        rhs: Vector<N, Sem, ScaleRhs, Space>,
    ) where
        Sem: FreelyMutable,
        ScaleRhs: ScaleSeal,
        Space: SpaceSeal,
    {
        for i in 0..N {
            self[i] += rhs[i];
        }
    }

    // SubAssign: both vectors must be in the same space
    #[impl_ops_assign(SubAssign)]
    fn sub_assign<const N: usize, Sem, ScaleRhs, Space>(
        self: &mut Vector<N, Sem, _NonUnit, Space>,
        rhs: Vector<N, Sem, ScaleRhs, Space>,
    ) where
        Sem: FreelyMutable,
        ScaleRhs: ScaleSeal,
        Space: SpaceSeal,
    {
        for i in 0..N {
            self[i] -= rhs[i];
        }
    }

    // MulAssign: component-wise, both vectors must be in the same space
    #[impl_ops_assign(MulAssign)]
    fn mul_assign<const N: usize, Sem, ScaleRhs, Space>(
        self: &mut Vector<N, Sem, _NonUnit, Space>,
        rhs: Vector<N, Sem, ScaleRhs, Space>,
    ) where
        Sem: FreelyMutable,
        ScaleRhs: ScaleSeal,
        Space: SpaceSeal,
    {
        for i in 0..N {
            self[i] *= rhs[i];
        }
    }

    // DivAssign: component-wise, both vectors must be in the same space
    #[impl_ops_assign(DivAssign)]
    fn div_assign<const N: usize, Sem, ScaleRhs, Space>(
        self: &mut Vector<N, Sem, _NonUnit, Space>,
        rhs: Vector<N, Sem, ScaleRhs, Space>,
    ) where
        Sem: FreelyMutable,
        ScaleRhs: ScaleSeal,
        Space: SpaceSeal,
    {
        for i in 0..N {
            self[i] /= rhs[i];
        }
    }
}

mod impl_rand {
    use super::*;

    use rand::Rng;
    use rand::distr::{Distribution, StandardUniform};

    impl<const N: usize, Sem> Distribution<Vector<N, Sem, _NonUnit, _WorldSpace>> for StandardUniform {
        fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vector<N, Sem, _NonUnit, _WorldSpace> {
            Vector {
                components: rng.random(),
                phantom: PhantomData,
            }
        }
    }
}

mod impl_default {
    use super::*;

    impl<const N: usize, Sem, Scale: ScaleSeal, Space: SpaceSeal> Default
        for Vector<N, Sem, Scale, Space>
    {
        fn default() -> Self {
            Self {
                components: [0.0; N],
                phantom: PhantomData,
            }
        }
    }
}

mod impl_index {
    use std::ops::{
        Index, IndexMut, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
    };

    use super::*;

    impl<const N: usize, Sem, Scale, Space> IndexMut<usize> for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, index: usize) -> &mut Self::Output {
            &mut self.components[index]
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<usize> for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = Float;

        fn index(&self, index: usize) -> &Self::Output {
            &self.components[index]
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<Range<usize>> for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = [Float];

        fn index(&self, range: Range<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<RangeFrom<usize>> for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = [Float];

        fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<RangeTo<usize>> for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = [Float];

        fn index(&self, range: RangeTo<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<RangeFull> for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = [Float; N];

        fn index(&self, _range: RangeFull) -> &Self::Output {
            &self.components
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<RangeInclusive<usize>>
        for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = [Float];

        fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> Index<RangeToInclusive<usize>>
        for Vector<N, Sem, Scale, Space>
    where
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        type Output = [Float];

        fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
            &self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> IndexMut<Range<usize>> for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, range: Range<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> IndexMut<RangeFrom<usize>> for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, range: RangeFrom<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> IndexMut<RangeTo<usize>> for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, range: RangeTo<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> IndexMut<RangeFull> for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, _range: RangeFull) -> &mut Self::Output {
            &mut self.components
        }
    }

    impl<const N: usize, Sem, Scale, Space> IndexMut<RangeInclusive<usize>>
        for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, range: RangeInclusive<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }

    impl<const N: usize, Sem, Scale, Space> IndexMut<RangeToInclusive<usize>>
        for Vector<N, Sem, Scale, Space>
    where
        Sem: FreelyMutable,
        Scale: ScaleSeal,
        Space: SpaceSeal,
    {
        fn index_mut(&mut self, range: RangeToInclusive<usize>) -> &mut Self::Output {
            &mut self.components[range]
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::Direction3;

    #[test]
    fn debug_test() {
        let v = Direction3::new(1.0, 2.0, 3.0);

        dbg!(v);
    }
}
