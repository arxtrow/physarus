/// Extracts and rearranges vector components by name.
///
/// # Examples
/// ```
/// use physarus::prelude::*;
///
/// let v = Direction3::new(1.0, 2.0, 3.0);
/// let xy = swizzle!(v => x, y);      // Vec2[1.0, 2.0]
/// let zyx = swizzle!(v => 2, 0, x);  // Vec3[3.0, 1.0, 1.0]
/// let zxy = swizzle!(v => z, x, 1); // Vector::<6, _Direction, _NonUnit>[3.0, 1.0, 2.0]
/// let xx = swizzle!(v => x, x);      // Vec2[1.0, 1.0]
/// let xyzzyx = swizzle!(v => x, y, z, z, y, x); // Vector::<6, _Direction, _NonUnit>[1.0, 2.0, 3.0, 3.0, 2.0, 1.0]
/// ```
#[macro_export]
macro_rules! swizzle {
    ($vec:expr => $($comp:tt),+ $(,)?) => {{
        #[allow(non_upper_case_globals)]
        {
            const x: usize = 0;
            const y: usize = 1;
            const z: usize = 2;
            const w: usize = 3;

            Direction::from_array([$($vec[$comp]),+])
        }
    }};
}
