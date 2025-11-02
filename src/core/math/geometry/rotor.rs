use std::marker::PhantomData;

use crate::core::{Float, Rad, Unit3, Vec3, Vec4, vec_states::*};

pub type Rotor<Space: SpaceSeal = _WorldSpace> = Vec4<_Rotor, _Unit, Space>;

impl<Space: SpaceSeal> Rotor<Space> {
    pub fn r0(&self) -> Float {
        self[0]
    }

    pub fn r12(&self) -> Float {
        self[1]
    }

    pub fn r23(&self) -> Float {
        self[2]
    }

    pub fn r31(&self) -> Float {
        self[3]
    }

    /// Params:
    /// `u` and `w`: non-parallel vectors defining the rotation plane and angle of rotation.
    ///
    /// Returns: normalized rotor encoding the rotation by the angle between `u` and `w` in plane of rotation defined by `u` and `w`.
    pub fn from_vecs(u: &Unit3<Space>, w: &Unit3<Space>) -> Self {
        debug_assert!(!u.is_parallel(w));

        let half = (u + w).normalize();

        let cos_half = u.dot(&half);
        let wedge = u.wedge(&half);

        let (r12, r23, r31) = (wedge.x(), wedge.y(), wedge.z());

        let rotor = Vec4 {
            components: [cos_half, r12, r23, r31],
            phantom: PhantomData,
        };

        rotor.normalize()
    }

    /// Params:
    /// `u`, `w`: non-parallel vectors defining the rotation plane
    /// `theta`: Rotation angle in radians
    ///
    /// Returns: normalized rotor encoding rotation by `theta` in the plane of `u` and `w`
    pub fn from_vecs_angle(u: &Unit3<Space>, w: &Unit3<Space>, theta: Rad) -> Self {
        debug_assert!(!u.is_parallel(w));

        let (sin_half, cos_half) = (theta / 2.0).sin_cos();
        let b = u.wedge(w);

        let r12 = sin_half * b.x();
        let r23 = sin_half * b.y();
        let r31 = sin_half * b.z();

        let rotor = Vec4 {
            components: [cos_half, r12, r23, r31],
            phantom: PhantomData,
        };

        rotor.normalize()
    }

    /// Params:
    /// `axis`: Unit vector defining rotation axis
    /// `theta`: Rotation angle in radians
    ///
    /// Returns: Rotor encoding rotation by `theta` around `axis`
    ///
    /// Formula: R = cos(θ/2) + sin(θ/2) * (n₃e₁₂ + n₂e₃₁ + n₁e₂₃ )
    pub fn from_axis_angle(axis: &Unit3<Space>, theta: Rad) -> Self {
        let (sin_half, cos_half) = (theta / 2.0).sin_cos();

        // The bivector components encode the plane perpendicular to the axis
        // In 3D GA: axis (n₁, n₂, n₃) is dual to bivector (n₃e₁₂, n₁e₂₃, n₂e₃₁)
        let r12 = sin_half * axis.z();
        let r23 = sin_half * axis.x();
        let r31 = sin_half * axis.y();

        let rotor = Vec4 {
            components: [cos_half, r12, r23, r31],
            phantom: PhantomData,
        };

        rotor.normalize()
    }

    /// Applies rotor rotation to vector using geometric algebra sandwich product (R v R~).
    /// Rotates input by the angle and plane encoded in the rotor, preserving vector magnitude.
    /// Params: `v: Vec3<SemRhs, ScaleRhs>` — vector to rotate.
    /// Returns rotated vector with same semantic and scale.
    pub fn rotate<SemRhs: GeoSeal, ScaleRhs: ScaleSeal>(
        &self,
        v: &Vec3<SemRhs, ScaleRhs, Space>,
    ) -> Vec3<SemRhs, ScaleRhs, Space> {
        let [r0, r12, r23, r31] = &self[..];
        let [vx, vy, vz] = &v[..];

        // First product: S = R v
        let (s1, s2, s3, s123) = (
            r0 * vx - r12 * vy + r31 * vz,
            r0 * vy - r23 * vz + r12 * vx,
            r0 * vz - r31 * vx + r23 * vy,
            r23 * vx + r31 * vy + r12 * vz,
        );

        // Second product: S R~
        let xyz = [
            s1 * r0 - s2 * r12 + s3 * r31 + s123 * r23,
            s2 * r0 + s1 * r12 - s3 * r23 + s123 * r31,
            s3 * r0 + s2 * r23 - s1 * r31 + s123 * r12,
        ];

        // SAFETY<Scale>: rotation does not change scale
        unsafe { Vec3::raw(xyz) }
    }

    /// Combines two rotors by geometric product
    /// Params: `self` - left (post) rotor; `rhs` - right (pre) rotor, applied first.
    /// Returns: new rotor representing the net rotation by `rhs` and then `self`.
    #[rustfmt::skip]
    pub fn combine(&self, rhs: &Self) -> Self {
        let scalar =
            self.r0() * rhs.r0()
            - self.r12() * rhs.r12()
            - self.r23() * rhs.r23()
            - self.r31() * rhs.r31();

        let xy =
            self.r0() * rhs.r12()
            + self.r12() * rhs.r0()
            - self.r23() * rhs.r31()
            + self.r31() * rhs.r23();

        let yz =
            self.r0() * rhs.r23()
            + self.r12() * rhs.r31()
            + self.r23() * rhs.r0()
            - self.r31() * rhs.r12();

        let zx =
            self.r0() * rhs.r31()
            - self.r12() * rhs.r23()
            + self.r23() * rhs.r12()
            + self.r31() * rhs.r0();

        let rotor = Vec4 {
            components: [scalar, xy, yz, zx],
            phantom: PhantomData,
        };

        rotor.normalize()
    }

    pub fn reverse(&self) -> Self {
        // SAFETY<Unit>: negation does not change length
        unsafe { Self::raw([self.r0(), -self.r12(), -self.r23(), -self.r31()]) }
    }
}

#[cfg(test)]
mod rotor_tests {
    use crate::core::{Direction3, PI, Rad, Vec3, traits::*};

    use super::*;

    #[test]
    fn test_rotation_smoke_test() {
        let sqrt2_inv = 1.0 / 2.0_f32.sqrt();

        let tests = [
            // Simple 90° rotations
            (Vec3::Z_AXIS, PI / 2.0, Vec3::X_AXIS, [0.0, 1.0, 0.0]),
            (Vec3::X_AXIS, PI / 2.0, Vec3::Y_AXIS, [0.0, 0.0, 1.0]),
            (Vec3::Y_AXIS, PI / 2.0, Vec3::Z_AXIS, [1.0, 0.0, 0.0]),
            // 180° and negative angles
            (Vec3::Z_AXIS, PI, Vec3::X_AXIS, [-1.0, 0.0, 0.0]),
            (Vec3::Z_AXIS, -PI / 2.0, Vec3::X_AXIS, [0.0, -1.0, 0.0]),
            (Vec3::Z_AXIS, 2.0 * PI, Vec3::X_AXIS, [1.0, 0.0, 0.0]),
            // 45° rotations (non-mainstream)
            (
                Vec3::Z_AXIS,
                PI / 4.0,
                Vec3::X_AXIS,
                [sqrt2_inv, sqrt2_inv, 0.0],
            ),
            (
                Vec3::X_AXIS,
                PI / 4.0,
                Vec3::Y_AXIS,
                [0.0, sqrt2_inv, sqrt2_inv],
            ),
            // 30° rotation around X-axis applied to Y-axis
            (Vec3::X_AXIS, PI / 6.0, Vec3::Y_AXIS, [0.0, 0.866025, 0.5]),
            // 60° rotation around Z-axis applied to X-axis
            (Vec3::Z_AXIS, PI / 3.0, Vec3::X_AXIS, [0.5, 0.866025, 0.0]),
            // Arbitrary axis: (1,1,1) normalized, 120° rotation on (1,0,0)
            (
                Direction3::new(1.0, 1.0, 1.0).normalize(),
                2.0 * PI / 3.0,
                Vec3::X_AXIS,
                [0.0, 1.0, 0.0],
            ),
            // Arbitrary axis: (1,1,0) normalized, 180° on (1,0,0)
            (
                Direction3::new(1.0, 1.0, 0.0).normalize(),
                PI,
                Vec3::X_AXIS,
                [0.0, 1.0, 0.0],
            ),
            // 72° rotation (pentagon angle) around Z on X
            (
                Vec3::Z_AXIS,
                PI * 2.0 / 5.0,
                Vec3::X_AXIS,
                [0.309017, 0.951057, 0.0],
            ),
            // Small angle: 0.1 radians (~5.7°) around Y-axis on X-axis
            (Vec3::Y_AXIS, 0.1, Vec3::X_AXIS, [0.995004, 0.0, -0.099833]),
            // Oblique axis (2,1,2) normalized, 90° on (1,0,0)
            (
                Direction3::new(2.0, 1.0, 2.0).normalize(),
                PI / 2.0,
                Vec3::X_AXIS,
                [0.444444, 0.888889, 0.111111],
            ),
        ];

        for (i, (axis, angle, input, expected)) in tests.iter().enumerate() {
            let rotor = Rotor::from_axis_angle(axis, Rad::new(*angle));
            let result = rotor.rotate(input);

            assert!(
                result.x().eq_eps(&expected[0])
                    && result.y().eq_eps(&expected[1])
                    && result.z().eq_eps(&expected[2]),
                "Test {} failed: expected [{:.6}, {:.6}, {:.6}], got [{:.6}, {:.6}, {:.6}]",
                i,
                expected[0],
                expected[1],
                expected[2],
                result.x(),
                result.y(),
                result.z()
            );
        }
    }

    #[test]
    fn test_axis_angle_rotation() {
        // Rotate x-axis by 90° around z-axis → should point along y-axis
        let z_axis = Vec3::Z_AXIS;
        let rotor = Rotor::from_axis_angle(&z_axis, Rad::new(PI / 2.0));

        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor.rotate(&x_vec);

        assert!(result.x().is_zero_eps(), "X component should be ~0");
        assert!(result.y().eq_eps(&1.0), "Y component should be 1");
        assert!(result.z().is_zero_eps(), "Z component should be ~0");
    }

    #[test]
    fn test_vecs_angle_rotation() {
        // Rotate in xy-plane from x-axis toward y-axis by 90°
        let x_axis = Vec3::X_AXIS;
        let y_axis = Vec3::Y_AXIS;
        let rotor = Rotor::from_vecs_angle(&x_axis, &y_axis, Rad::new(PI / 2.0));

        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor.rotate(&x_vec);

        assert!(result.x().is_zero_eps(), "Should rotate to y-axis: x ≈ 0");
        assert!(result.y().eq_eps(&1.0), "Should rotate to y-axis: y ≈ 1");
        assert!(result.z().is_zero_eps(), "Should rotate to y-axis: z ≈ 0");
    }

    #[test]
    fn test_from_vecs_full_rotation() {
        // Rotate from x-axis to y-axis (90° implicitly)
        let x_axis = Vec3::X_AXIS;
        let y_axis = Vec3::Y_AXIS;
        let rotor = Rotor::from_vecs(&x_axis, &y_axis);

        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor.rotate(&x_vec);

        assert!(result.x().is_zero_eps(), "Should map x to y");
        assert!(result.y().eq_eps(&1.0), "Should map x to y");
        assert!(result.z().is_zero_eps(), "Should map x to y");
    }

    #[test]
    fn test_three_methods_equivalent() {
        // All three methods should produce the same rotation
        let x_axis = Vec3::X_AXIS;
        let y_axis = Vec3::Y_AXIS;
        let z_axis = Vec3::Z_AXIS;

        let rotor1 = Rotor::from_axis_angle(&z_axis, Rad::new(PI / 2.0));
        let rotor2 = Rotor::from_vecs_angle(&x_axis, &y_axis, Rad::new(PI / 2.0));
        let rotor3 = Rotor::from_vecs(&x_axis, &y_axis);

        let test_vec = Direction3::new(1.0, 0.0, 0.0);

        let result1 = rotor1.rotate(&test_vec);
        let result2 = rotor2.rotate(&test_vec);
        let result3 = rotor3.rotate(&test_vec);

        // Compare component-wise for better error messages
        assert!(result1.x().eq_eps(&result2.x()), "Method 1 vs 2: X differs");
        assert!(result1.y().eq_eps(&result2.y()), "Method 1 vs 2: Y differs");
        assert!(result1.z().eq_eps(&result2.z()), "Method 1 vs 2: Z differs");

        assert!(result1.x().eq_eps(&result3.x()), "Method 1 vs 3: X differs");
        assert!(result1.y().eq_eps(&result3.y()), "Method 1 vs 3: Y differs");
        assert!(result1.z().eq_eps(&result3.z()), "Method 1 vs 3: Z differs");
    }

    #[test]
    fn test_45_degree_rotation() {
        // Rotate by 45° in xy-plane
        let x_axis = Vec3::X_AXIS;
        let y_axis = Vec3::Y_AXIS;
        let rotor = Rotor::from_vecs_angle(&x_axis, &y_axis, Rad::new(PI / 4.0));

        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor.rotate(&x_vec);

        let expected = 1.0 / 2.0_f32.sqrt();
        assert!(result.x().eq_eps(&expected), "X should be 1/√2");
        assert!(result.y().eq_eps(&expected), "Y should be 1/√2");
        assert!(result.z().is_zero_eps(), "Z should be 0");
    }

    #[test]
    fn test_identity_rotation() {
        // Rotating by 0° should do nothing
        let x_axis = Vec3::X_AXIS;
        let rotor = Rotor::from_axis_angle(&x_axis, Rad::new(0.0));

        let test_vec = Direction3::new(1.0, 2.0, 3.0).normalize();
        let result = rotor.rotate(&test_vec);

        assert!(
            (result.x() - test_vec.x()).is_zero_eps(),
            "Identity rotation should preserve X"
        );
        assert!(
            (result.y() - test_vec.y()).is_zero_eps(),
            "Identity rotation should preserve Y"
        );
        assert!(
            (result.z() - test_vec.z()).is_zero_eps(),
            "Identity rotation should preserve Z"
        );
    }

    #[test]
    fn test_180_degree_rotation() {
        // Rotate by 180° around z-axis: (1,0,0) → (-1,0,0)
        let z_axis = Vec3::Z_AXIS;
        let rotor = Rotor::from_axis_angle(&z_axis, Rad::new(PI));

        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor.rotate(&x_vec);

        assert!(result.x().eq_eps(&-1.0), "X should be -1");
        assert!(result.y().is_zero_eps(), "Y should be 0");
        assert!(result.z().is_zero_eps(), "Z should be 0");
    }

    #[test]
    fn test_rotor_composition() {
        // Compose two 45° rotations to get 90°
        let z_axis = Vec3::Z_AXIS;
        let rotor45 = Rotor::from_axis_angle(&z_axis, Rad::new(PI / 4.0));
        let rotor_combined = rotor45.combine(&rotor45);

        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor_combined.rotate(&x_vec);

        assert!(result.x().is_zero_eps(), "Combined rotation: X → 0");
        assert!(result.y().eq_eps(&1.0), "Combined rotation: Y → 1");
        assert!(result.z().is_zero_eps(), "Combined rotation: Z → 0");
    }

    #[test]
    fn test_reverse_rotor() {
        // Rotor and its reverse should undo each other
        let z_axis = Vec3::Z_AXIS;
        let rotor = Rotor::from_axis_angle(&z_axis, Rad::new(PI / 3.0));
        let rotor_rev = rotor.reverse();

        let test_vec = Direction3::new(1.0, 2.0, 3.0).normalize();
        let rotated = rotor.rotate(&test_vec);
        let back = rotor_rev.rotate(&rotated);

        assert!(
            (back.x() - test_vec.x()).is_zero_eps(),
            "Reverse should undo rotation: X component"
        );
        assert!(
            (back.y() - test_vec.y()).is_zero_eps(),
            "Reverse should undo rotation: Y component"
        );
        assert!(
            (back.z() - test_vec.z()).is_zero_eps(),
            "Reverse should undo rotation: Z component"
        );
    }

    #[test]
    fn test_arbitrary_axis_rotation() {
        // Rotate around arbitrary axis (1,1,1)
        let axis = Vec3::new(1.0, 1.0, 1.0).normalize();
        let rotor = Rotor::from_axis_angle(&axis, Rad::new(2.0 * PI / 3.0));

        // Rotating (1,0,0) by 120° around (1,1,1) should give (0,1,0)
        let x_vec = Direction3::new(1.0, 0.0, 0.0);
        let result = rotor.rotate(&x_vec);

        assert!(result.x().is_zero_eps(), "120° rotation cycles x→y: X ≈ 0");
        assert!(result.y().eq_eps(&1.0), "120° rotation cycles x→y: Y ≈ 1");
        assert!(result.z().is_zero_eps(), "120° rotation cycles x→y: Z ≈ 0");
    }

    #[test]
    fn test_preserves_magnitude() {
        // Rotation should preserve vector magnitude
        let axis = Direction3::new(1.0, 2.0, 3.0).normalize();
        let rotor = Rotor::from_axis_angle(&axis, Rad::new(0.8));

        let test_vec = Direction3::new(3.0, 4.0, 5.0);
        let original_len = test_vec.len();
        let result = rotor.rotate(&test_vec);
        let result_len = result.len();

        assert!(
            original_len.eq_eps(&result_len),
            "Rotation should preserve magnitude: {} vs {}",
            original_len,
            result_len
        );
    }

    #[test]
    fn test_orthogonal_vectors_remain_orthogonal() {
        // Rotating orthogonal vectors should keep them orthogonal
        let z_axis = Unit3::Z_AXIS;
        let rotor = Rotor::from_axis_angle(&z_axis, Rad::new(0.7));

        let v1 = Vec3::X_AXIS;
        let v2 = Vec3::Y_AXIS;

        // Original dot product should be 0 (orthogonal)
        assert!(
            v1.dot(&v2).is_zero_eps(),
            "Original vectors should be orthogonal"
        );

        let r1 = rotor.rotate(&v1);
        let r2 = rotor.rotate(&v2);

        // Rotated vectors should still be orthogonal
        assert!(
            r1.dot(&r2).is_zero_eps(),
            "Rotated vectors should remain orthogonal"
        );
    }
}
