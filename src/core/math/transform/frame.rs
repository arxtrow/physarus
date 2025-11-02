use crate::prelude::*;

/// Frame represented by 3 orthonormal vectors
/// Righthanded: x.cross(y) = z
/// Given 3 orthonormal vectors x, y, z, the matrix F that transforms vectors INTO their space is
/// F = [
///     --- x ---
///     --- y ---
///     --- z ---
/// ]
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Frame {
    x: Unit3,
    y: Unit3,
    z: Unit3,
}

impl Default for Frame {
    fn default() -> Self {
        Self {
            x: Unit3::X_AXIS,
            y: Unit3::Y_AXIS,
            z: Unit3::Z_AXIS,
        }
    }
}

impl Frame {
    pub fn new(x: Unit3, y: Unit3, z: Unit3) -> Self {
        debug_assert!(x.is_perpendicular(&y) && y.is_perpendicular(&z) && z.is_perpendicular(&x));
        debug_assert!(x.cross(&y).dot(&z) > 0.0, "Frame must be right-handed");

        Self { x, y, z }
    }

    pub fn from_xz(x: Unit3, z: Unit3) -> Self {
        let y = z.cross(&x).normalize();

        Self { x, y, z }
    }

    pub fn from_xy(x: Unit3, y: Unit3) -> Self {
        let z = x.cross(&y).normalize();

        Self { x, y, z }
    }

    pub fn from_z(z: Unit3) -> Self {
        let (x, y, z) = z.coordinate_system();

        Self { x, y, z }
    }

    pub fn from_y(y: Unit3) -> Self {
        let (z, x, y) = y.coordinate_system();

        Self { x, y, z }
    }

    pub fn from_x(x: Unit3) -> Self {
        let (y, z, x) = x.coordinate_system();

        Self { x, y, z }
    }

    pub fn from_normal(n: &Normal3) -> Self {
        let (x, y, z) = n.coordinate_system();

        Self { x, y, z }
    }

    /// works with `Normal3` too: F is orthonormal, so `n' = (F⁻¹)ᵀ n = (Fᵀ)ᵀ n = F n` and det(F) = 1
    pub fn to_local(
        &self,
        vec: &Vec3<impl GeoSeal, impl ScaleSeal, _WorldSpace>,
    ) -> Direction3<_FrameSpace> {
        Direction3::new(self.x.dot(&vec), self.y.dot(&vec), self.z.dot(&vec)).in_space()
    }

    /// `F⁻¹ = Fᵀ = [columns: x y z]`
    /// Represent `(F⁻¹ vec)` as linear combination of columns with coefficients from vec
    /// Works with `Normal3` since `F` is orthonormal `(det(F) = 1)`
    pub fn from_local(
        &self,
        vec: &Vec3<impl GeoSeal, impl ScaleSeal, _FrameSpace>,
    ) -> Direction3<_WorldSpace> {
        self.x * vec.x() + self.y * vec.y() + self.z * vec.z()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_constructors_are_right_handed() {
        let x = Unit3::X_AXIS;
        let y = Unit3::Y_AXIS;
        let z = Unit3::Z_AXIS;

        let frame_default = Frame::default();
        assert!(
            frame_default
                .x
                .cross(&frame_default.y)
                .dot(&frame_default.z)
                > 0.0
        );

        let frame_xy = Frame::from_xy(x, y);
        assert!(frame_xy.x.cross(&frame_xy.y).dot(&frame_xy.z) > 0.0);

        let frame_xz = Frame::from_xz(x, z);
        assert!(frame_xz.x.cross(&frame_xz.y).dot(&frame_xz.z) > 0.0);

        let frame_z = Frame::from_z(z);
        assert!(frame_z.x.cross(&frame_z.y).dot(&frame_z.z) > 0.0);
    }

    #[test]
    fn frame_with_arbitrary_up_vector() {
        let z = rand::random::<Direction3>().normalize();
        let frame = Frame::from_z(z);

        // Should be orthonormal and right-handed
        assert!(frame.x.is_perpendicular(&frame.y));
        assert!(frame.y.is_perpendicular(&frame.z));
        assert!(frame.z.is_perpendicular(&frame.x));
        assert!(frame.x.cross(&frame.y).dot(&frame.z) > 0.0);
    }

    #[test]
    fn frame_back_and_forth_preserve() {
        let dir = rand::random::<Direction3>();
        let z = rand::random::<Direction3>().normalize();
        let frame = Frame::from_z(z);

        assert!(dir.eq_eps(&frame.from_local(&frame.to_local(&dir))))
    }

    #[test]
    fn to_local_preserve_dot() {
        let dir1 = rand::random::<Direction3>();
        let dir2 = rand::random::<Direction3>();
        let z = rand::random::<Direction3>().normalize();
        let frame = Frame::from_z(z);

        let dir1_local = frame.to_local(&dir1);
        let dir2_local = frame.to_local(&dir2);

        assert!((dir1.dot(&dir2)).eq_eps(&dir1_local.dot(&dir2_local)))
    }

    #[test]
    fn to_local_works_with_normals() {
        let v1_world = rand::random::<Direction3>();
        let v2_world = rand::random::<Direction3>();
        let n_world = Normal3::from_cross(&v1_world, &v2_world);

        let frame = Frame::from_z(rand::random::<Direction3>().normalize());

        assert!(
            frame
                .to_local(&n_world)
                .is_perpendicular(&frame.to_local(&v1_world))
        );

        assert!(
            frame
                .to_local(&n_world)
                .is_perpendicular(&frame.to_local(&v2_world))
        );
    }

    #[test]
    fn from_local_works_with_normals() {
        let v1_local = rand::random::<Direction3>().in_space::<_FrameSpace>();
        let v2_local = rand::random::<Direction3>().in_space::<_FrameSpace>();
        let n_local = Normal3::from_cross(&v1_local, &v2_local).in_space::<_FrameSpace>();

        let frame = Frame::from_z(rand::random::<Direction3>().normalize());

        assert!(
            frame
                .from_local(&n_local)
                .is_perpendicular(&frame.from_local(&v1_local))
        );

        assert!(
            frame
                .from_local(&n_local)
                .is_perpendicular(&frame.from_local(&v2_local))
        );
    }

    #[test]
    fn from_xyz_all_preserve_input_axis() {
        let x_input = Unit3::X_AXIS;
        let y_input = Unit3::Y_AXIS;
        let z_input = Unit3::Z_AXIS;

        let frame_x = Frame::from_x(x_input);
        let frame_y = Frame::from_y(y_input);
        let frame_z = Frame::from_z(z_input);

        // Check the correct axis is preserved
        assert_eq!(frame_x.x, x_input);
        assert_eq!(frame_y.y, y_input);
        assert_eq!(frame_z.z, z_input);

        // Check all are right-handed
        assert!(frame_x.x.cross(&frame_x.y).dot(&frame_x.z) > 0.0);
        assert!(frame_y.x.cross(&frame_y.y).dot(&frame_y.z) > 0.0);
        assert!(frame_z.x.cross(&frame_z.y).dot(&frame_z.z) > 0.0);
    }
}
