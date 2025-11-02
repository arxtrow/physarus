use std::{marker::PhantomData, ops::Mul, simd::Simd};

use opimps::impl_ops;

use crate::prelude::*;

/// A transformation matrix with tracked source and destination coordinate spaces.
/// Generic over `From` and `To` spaces types to ensure compile-time correctness when composing
/// transformations between different coordinate systems.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Transform<From: SpaceSeal = _WorldSpace, To: SpaceSeal = _WorldSpace> {
    m: Mat4,
    m_inv: Mat4,

    phantom: PhantomData<(From, To)>,
}

impl Transform<_WorldSpace, _WorldSpace> {
    /// Create a world-to-world transform (most common case).
    /// ***The caller must ensure that `m_inv` is the correct matrix inverse of `m`.***
    pub fn new(m: Mat4, m_inv: Mat4) -> Self {
        Self {
            m,
            m_inv,
            phantom: PhantomData,
        }
    }

    pub fn identity() -> Self {
        let m = Mat4::identity();
        let m_inv = Mat4::identity();

        Transform::new_from_to(m, m_inv)
    }

    pub fn from_translation(v: &Direction3) -> Self {
        let m = mat4!(
            1.0, 0.0, 0.0, v.x(),
            0.0, 1.0, 0.0, v.y(),
            0.0, 0.0, 1.0, v.z(),
            0.0, 0.0, 0.0, 1.0,
        );

        let m_inv = mat4!(
            1.0, 0.0, 0.0, -v.x(),
            0.0, 1.0, 0.0, -v.y(),
            0.0, 0.0, 1.0, -v.z(),
            0.0, 0.0, 0.0, 1.0,
        );

        Self::new(m, m_inv)
    }

    pub fn from_scale(x: Float, y: Float, z: Float) -> Self {
        let m = mat4!(
            x,   0.0, 0.0, 0.0,
            0.0, y,   0.0, 0.0,
            0.0, 0.0, z,   0.0,
            0.0, 0.0, 0.0, 1.0,
        );

        let m_inv = mat4!(
            1.0 / x, 0.0,     0.0,     0.0,
            0.0,     1.0 / y, 0.0,     0.0,
            0.0,     0.0,     1.0 / z, 0.0,
            0.0,     0.0,     0.0,     1.0,
        );

        Self::new(m, m_inv)
    }

    pub fn from_rotor(rotor: &Rotor) -> Self {
        let x_axis = rotor.rotate(&Unit3::X_AXIS);
        let y_axis = rotor.rotate(&Unit3::Y_AXIS);
        let z_axis = rotor.rotate(&Unit3::Z_AXIS);

        let m = mat4!(
            x_axis.x(), y_axis.x(), z_axis.x(), 0.0,
            x_axis.y(), y_axis.y(), z_axis.y(), 0.0,
            x_axis.z(), y_axis.z(), z_axis.z(), 0.0,
            0.0,        0.0,        0.0,        1.0,
        );

        let m_inv = m.transpose();

        Self::new(m, m_inv)
    }
    //
    // pub fn around_x(theta: Rad) -> Self {
    //     let (sin, cos) = theta.sin_cos();
    //
    //     let m = mat4!(
    //     todo!()
    // );
    //
    //     let m_inv = m.transpose();
    //
    //     Self::new(m, m_inv)
    // }
    //
    // pub fn around_y(theta: Rad) -> Self {
    //     let (sin, cos) = theta.sin_cos();
    //
    //     let m = mat4!(
    //     todo!()
    // );
    //
    //     let m_inv = m.transpose();
    //
    //     Self::new(m, m_inv)
    // }
    //
    // pub fn around_z(theta: Rad) -> Self {
    //     let (sin, cos) = theta.sin_cos();
    //
    //     let m = mat4!(
    //     todo!()
    // );
    //
    //     let m_inv = m.transpose();
    //
    //     Self::new(m, m_inv)
    // }

    pub fn has_scale(&self) -> bool {
        let x_len_sqr = self.transform(&Unit3::X_AXIS).len_sqr();
        let y_len_sqr = self.transform(&Unit3::Y_AXIS).len_sqr();
        let z_len_sqr = self.transform(&Unit3::Z_AXIS).len_sqr();

        !x_len_sqr.eq_eps(&1.0) || !y_len_sqr.eq_eps(&1.0) || !x_len_sqr.eq_eps(&1.0)
    }

    pub fn swaps_handedness(&self) -> bool {
        self.m.upper3x3().det() < 0.0
    }
}

impl<From: SpaceSeal, To: SpaceSeal> Transform<From, To> {
    /// Create a transform between arbitrary coordinate spaces.
    ///
    /// Requires both forward (m) and inverse (m_inv) matrices.
    /// The caller must ensure m_inv is actually the inverse of m.
    pub fn new_from_to(m: Mat4, m_inv: Mat4) -> Self {
        Self {
            m,
            m_inv,
            phantom: PhantomData,
        }
    }

    pub fn transform<T: Transformable<From = From>>(&self, value: &T) -> T::Output<To> {
        value.apply_transform(self)
    }

    pub fn inverse_transform<T: Transformable<From = From>>(&self, value: &T) -> T::Output<From> {
        value.apply_inverse(self)
    }

    pub fn inverse(&self) -> Transform<To, From> {
        Transform::new_from_to(self.m_inv, self.m)
    }

    pub fn is_identity(&self) -> bool {
        self.m.is_identity()
    }

    pub fn matrix(&self) -> Mat4 {
        self.m.clone()
    }

    pub fn inverse_matrix(&self) -> Mat4 {
        self.m_inv.clone()
    }
}

impl Transform<_WorldSpace, _CameraSpace> {
    pub fn look_at(
        pos: &Point3<_WorldSpace>,
        target: &Point3<_WorldSpace>,
        up: Unit3<_WorldSpace>,
    ) -> Self {
        let z = (pos - target).normalize(); // camera looks toward -z, +z points to screen 
        let x = up.cross(&z).normalize(); // right
        let y = z.cross(&x).normalize(); // up

        // Camera-to-World transformation: [R | pos]
        let camera_to_world = mat4!(
            x.x(), y.x(), z.x(), pos.x(),
            x.y(), y.y(), z.y(), pos.y(),
            x.z(), y.z(), z.z(), pos.z(),
            0.0,   0.0,   0.0,    1.0,
        );

        // World-to-Camera (view) transformation = [ R^T | -R^T pos ]
        // -R^T pos: negative because we want camera at origin
        // and translation vector precomputed in camera space since we apply translation
        // to vector in camera space: [ R^T | -R^T*pos ] * v = (R^T*v) + (-R^T*pos))
        let world_to_camera = mat4!(
            x.x(), x.y(), x.z(), -x.dot(&pos.into_direction()),
            y.x(), y.y(), y.z(), -y.dot(&pos.into_direction()),
            z.x(), z.y(), z.z(), -z.dot(&pos.into_direction()),
            0.0,   0.0,   0.0,    1.0,
        );

        Transform::new_from_to(world_to_camera, camera_to_world)
    }
}

/// Compose two transforms: self ∘ rhs, producing from_rhs → to_self.
/// Inverses composed in reverse order to comply with (A B)⁻¹ = B⁻¹A⁻¹
#[impl_ops(Mul)]
fn mul<Mid, ToLhs, FromRhs>(
    self: Transform<Mid, ToLhs>,
    rhs: Transform<FromRhs, Mid>,
) -> Transform<FromRhs, ToLhs>
where
    Mid: SpaceSeal,
    ToLhs: SpaceSeal,
    FromRhs: SpaceSeal,
{
    return Transform::new_from_to(self.m * rhs.m, rhs.m_inv * self.m_inv);
}

pub trait Transformable {
    type From: SpaceSeal;

    type Output<To: SpaceSeal>;

    /// Applies Transform in semantic correct way
    fn apply_transform<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<To>;

    /// Applies Inverse Transform in semantic correct way
    fn apply_inverse<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<Self::From>;
}

impl<Scale: ScaleSeal, From: SpaceSeal> Transformable for Vec3<_Direction, Scale, From> {
    type From = From;
    type Output<To: SpaceSeal> = Direction3<To>;

    /// Transform `Direction` with `M₃ₓ₃`
    fn apply_transform<To: SpaceSeal>(&self, transform: &Transform<From, To>) -> Self::Output<To> {
        let m3x3 = transform.m.upper3x3();

        let tv = m3x3 * self;

        Direction3::from_array_in(tv.to_array())
    }

    /// Inverse Transform `Direction` with `(M₃ₓ₃)⁻¹`
    fn apply_inverse<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<Self::From> {
        let m_inv3x3 = transform.m_inv.upper3x3();

        let iv = m_inv3x3 * self;

        Direction3::from_array_in(iv.to_array())
    }
}

impl<From: SpaceSeal> Transformable for Normal3<From> {
    type From = From;
    type Output<To: SpaceSeal> = Normal3<To>;

    /// Transform `Normal3` with `((M₃ₓ₃)⁻¹)ᵀ`
    fn apply_transform<To: SpaceSeal>(&self, transform: &Transform<From, To>) -> Self::Output<To> {
        // Normals are Hodge duals of bivectors
        // bivectors represent oriented planes with area and transform extensively (area-like)
        // then normals must transform inversely to preserve geometric relationships like perpendicularity with vectors in plane.
        //
        // Derivation: For tangent vector t in the plane and any vector v, we have
        // n·t = 0 and n·v = k. After transformation by M, we want:
        //   n'·(Mt) = n·t = 0  and  n'·(Mv) = n·v = k
        //
        // In matrix form for all vectors: n'ᵀM = nᵀ
        //   n'ᵀ = nᵀM⁻¹
        //   n'ᵀ = ((M⁻¹)ᵀn)ᵀ
        //   n' = (M⁻¹)ᵀn
        //
        // When M⁻¹ is not precomputed, use the cofactor matrix:
        //   cof(M) = adj(M)ᵀ = det(M)·(M⁻¹)ᵀ
        // so n' = cof(M)·n, then normalize (cancels det(M) factor)
        let tn = transform.m_inv.upper3x3().transpose() * self;

        // Normalize: (M⁻¹)ᵀ doesn't preserve length under non-uniform scaling/shear
        Normal3::new_normalize(tn.x(), tn.y(), tn.z())
    }

    /// Inverse Transform `Normal3` with `(M₃ₓ₃)ᵀ`
    fn apply_inverse<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<Self::From> {
        let tn = transform.m.upper3x3().transpose() * self;

        Normal3::new_normalize(tn.x(), tn.y(), tn.z())
    }
}

impl<From: SpaceSeal> Transformable for Point3<From> {
    type From = From;
    type Output<To: SpaceSeal> = Point3<To>;

    /// Transform `Point3` using full `M`
    ///
    /// Applies perspective division after transformation.
    fn apply_transform<To: SpaceSeal>(&self, transform: &Transform<From, To>) -> Self::Output<To> {
        let v = Direction4::<From>::new(self[0], self[1], self[2], 1.0);

        let tv = transform.m * v;
        let [xp, yp, zp, wp] = tv.to_array();

        // Perspective division
        Point3::new_array([xp, yp, zp]) / wp
    }

    /// Inverse Transform `Point3` using full `M⁻¹`
    ///
    /// Applies perspective division after transformation.
    fn apply_inverse<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<Self::From> {
        let v = Direction4::<From>::new(self[0], self[1], self[2], 1.0);

        let tv = transform.m_inv * v;
        let [xp, yp, zp, wp] = tv.to_array();

        // Perspective division
        Point3::new_array([xp, yp, zp]) / wp
    }
}

impl<From: SpaceSeal> Transformable for Ray3<From> {
    type From = From;
    type Output<To: SpaceSeal> = Ray3<To>;

    fn apply_transform<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<To> {
        let o = self.o.apply_transform(transform);
        let d = self.d.apply_transform(transform).normalize();

        let t_max = self.t_max.get();

        let len_sqr = d.len_sqr();
        if len_sqr >= 0.0 {
            // FIXME
        }

        Ray3::new(o, d, self.time, t_max)
    }

    fn apply_inverse<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<Self::From> {
        let o = self.o.apply_inverse(transform);
        let d = self.d.apply_inverse(transform).normalize();

        let t_max = self.t_max.get();

        let len_sqr = d.len_sqr();
        if len_sqr >= 0.0 {
            // FIXME
        }

        Ray3::new(o, d, self.time, t_max)
    }
}

impl<From: SpaceSeal> Transformable for Bounds3<From> {
    type From = From;
    type Output<To: SpaceSeal> = Bounds3<To>;

    fn apply_transform<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<To> {
        let mut bb = Bounds3::default();

        for i in 0..8 {
            let t_corner = self.corner(i).apply_transform(&transform);
            bb = bb.union_point(&t_corner);
        }

        bb
    }

    fn apply_inverse<To: SpaceSeal>(
        &self,
        transform: &Transform<Self::From, To>,
    ) -> Self::Output<Self::From> {
        todo!()
    }
}
