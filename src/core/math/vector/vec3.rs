use std::marker::PhantomData;

use crate::prelude::*;

pub type Vec3<Sem, Scale: ScaleSeal, Space: SpaceSeal> = Vector<3, Sem, Scale, Space>;

impl Vec3<_Direction, _NonUnit, _WorldSpace> {
    pub const ZERO: Self = Vec3 {
        components: [0.0, 0.0, 0.0],
        phantom: PhantomData,
    };

    pub const ONE: Self = Vec3 {
        components: [1.0, 1.0, 1.0],
        phantom: PhantomData,
    };
}

impl Unit3 {
    pub const I: Self = Self {
        components: [1.0, 0.0, 0.0],
        phantom: PhantomData,
    };

    pub const J: Self = Self {
        components: [0.0, 1.0, 0.0],
        phantom: PhantomData,
    };

    pub const K: Self = Self {
        components: [0.0, 0.0, 1.0],
        phantom: PhantomData,
    };

    pub const X_AXIS: Self = Self {
        components: [1.0, 0.0, 0.0],
        phantom: PhantomData,
    };

    pub const Y_AXIS: Self = Self {
        components: [0.0, 1.0, 0.0],
        phantom: PhantomData,
    };

    pub const Z_AXIS: Self = Self {
        components: [0.0, 0.0, 1.0],
        phantom: PhantomData,
    };
}

impl<Sem> Vec3<Sem, _NonUnit, _WorldSpace> {
    pub const fn new(x: Float, y: Float, z: Float) -> Self {
        Vector {
            components: [x, y, z],
            phantom: PhantomData,
        }
    }
}

impl<Sem, Scale, Space> Vec3<Sem, Scale, Space>
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

    pub fn z(&self) -> Float {
        self[2]
    }
}

impl<SemLhs, ScaleLhs, Space> Vec3<SemLhs, ScaleLhs, Space>
where
    SemLhs: GeoSeal,
    ScaleLhs: ScaleSeal,
    Space: SpaceSeal,
{
    pub fn cross(&self, rhs: &Vec3<impl GeoSeal, impl ScaleSeal, Space>) -> Direction3<Space> {
        Direction3::from_array_in([
            (self[1], rhs[2]).diff_of_mul((self[2], rhs[1])), // y*z - z*y
            (self[2], rhs[0]).diff_of_mul((self[0], rhs[2])), // z*x - x*z
            (self[0], rhs[1]).diff_of_mul((self[1], rhs[0])), // x*y - y*x
        ])
    }

    pub fn wedge(
        &self,
        rhs: &Vec3<impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Vec3<_Bivector, _NonUnit, Space> {
        let cross = self.cross(rhs);
        unsafe { Vec3::raw([cross.z(), cross.x(), cross.y()]) }
    }
}

impl<SemLhs, ScaleLhs, Space> Vec3<SemLhs, ScaleLhs, Space>
where
    SemLhs: GeoSeal,
    ScaleLhs: ScaleSeal,
    Space: SpaceSeal,
{
    pub fn angle(&self, rhs: &Vec3<impl GeoSeal, impl ScaleSeal, Space>) -> Rad {
        let cross = self.cross(rhs).len(); // == |self||rhs| sin(angle)
        let dot = self.dot(rhs); // == |self||rhs| cos(angle)

        Rad::new(cross.atan2(dot))
    }
}

pub type Direction3<Space: SpaceSeal = _WorldSpace> = Vec3<_Direction, _NonUnit, Space>;
pub type Unit3<Space: SpaceSeal = _WorldSpace> = Vec3<_Direction, _Unit, Space>;

impl<Space: SpaceSeal> Unit3<Space> {
    #[inline]
    pub fn into_normal3(&self) -> Normal3<Space> {
        unsafe { Normal3::raw(self.components) }
    }

    // https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf
    // with fixes from https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    //
    // Geometric Intuition:
    // Formulas taken from quaternion that rotates the canonical z-axis (0, 0, 1) to align with normal n
    // And same quaternion simultaneously rotates the x-axes to b1 and y-axes to b2.
    // This yields a complete change of basis {b1, b2, n} that is itself an orthonormal frame
    //
    /// Constructs right-handed orthonormal basis from normal using quaternion-rotation derived formulas.
    /// Returns: `(tangent1, tangent2, normal)` forming right-handed orthonormal coordinate frame
    pub fn coordinate_system(&self) -> (Unit3<Space>, Unit3<Space>, Unit3<Space>) {
        let z = *self;

        let s = (1.0 as Float).copysign(z.z());
        let a = -1.0 / (s + z.z());
        let b = z.x() * z.y() * a;

        // SAFETY<Unit>: b1 is unit and b2 is unit. Proved in frisvad paper
        let x = unsafe { Unit3::raw([1.0 + s * z.x() * z.x() * a, s * b, -s * z.x()]) };
        let y = unsafe { Unit3::raw([b, s + z.y() * z.y() * a, -z.y()]) };

        (x, y, z)
    }

    pub fn new_normalize(x: Float, y: Float, z: Float) -> Self {
        let len_sqr = x * x + y * y + z * z;

        debug_assert!(len_sqr > 0.0, "Cannot normalize zero-length vector");

        let len = len_sqr.sqrt();

        // SAFETY<Unit>: normalized by len
        unsafe { Self::raw([x / len, y / len, z / len]) }
    }

    pub fn from_cross(
        v1: &Vec3<impl GeoSeal, impl ScaleSeal, Space>,
        v2: &Vec3<impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Self {
        let cross = v1.cross(v2).normalize();

        // SAFETY<Unit>: normalized by len
        unsafe { Self::raw(cross.components) }
    }

    pub fn face_forward(&self, v: &Vec3<impl GeoSeal, impl ScaleSeal, Space>) -> Self {
        let u = *self;
        if self.dot(v) < 0.0 { -u } else { u }
    }
}

pub type Radiance3<Scale: ScaleSeal = _NonUnit, Space: SpaceSeal = _WorldSpace> =
    Vec3<_Radiance, Scale, Space>;

pub type Tangent3<Space: SpaceSeal = _WorldSpace> = Vec3<_Tangent, _Unit, Space>;
pub type Normal3<Space: SpaceSeal = _WorldSpace> = Vec3<_Normal, _Unit, Space>;

impl<Space: SpaceSeal> Normal3<Space> {
    #[inline]
    pub fn into_unit3(&self) -> Unit3<Space> {
        unsafe { Unit3::raw(self.components) }
    }

    pub fn new_normalize(x: Float, y: Float, z: Float) -> Self {
        Unit3::new_normalize(x, y, z).into_normal3()
    }

    pub fn from_cross(
        v1: &Vec3<impl GeoSeal, impl ScaleSeal, Space>,
        v2: &Vec3<impl GeoSeal, impl ScaleSeal, Space>,
    ) -> Self {
        Unit3::from_cross(v1, v2).into_normal3()
    }

    pub fn coordinate_system(&self) -> (Unit3<Space>, Unit3<Space>, Unit3<Space>) {
        self.into_unit3().coordinate_system()
    }

    pub fn face_forward(&self, v: &Vec3<impl GeoSeal, impl ScaleSeal, Space>) -> Self {
        self.into_unit3().face_forward(v).into_normal3()
    }
}

pub type Bitangent3<Space: SpaceSeal = _WorldSpace> = Vec3<_Bitangent, _Unit, Space>;

impl<Space: SpaceSeal> Bitangent3<Space> {
    pub fn from_normal_tangent(normal: &Normal3<Space>, tangent: &Tangent3<Space>) -> Self {
        debug_assert!(normal.is_perpendicular(&tangent));

        let cross = normal.cross(tangent).normalize();

        // SAFETY<Unit>: normalized
        unsafe { Vec3::raw(cross.components) }
    }
}

// FIXME How to deal with outside coords
// pub type Barycentric3 = Vec3<Barycentric, SumToOne>;
// impl Barycentric3 {
//     pub fn is_inside(&self) {
//         self.x() >= 0.0 && self.y() >= 0.0 && self.z() >= 0.0
//     }
// }
//
