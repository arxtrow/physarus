use crate::state_structs;

state_structs!(
    // Scale
    /// Any vector length (could be NaN, 1.0, 0.0, ...)
    _NonUnit,
    /// Vector.len() == 1.0
    _Unit,
    // Scale End

    // Semantic
    /// General-purpose geometric vector
    /// Serves as the default for `Vector` and the fallback result type when operations mix different semantic types (e.g., `Normal + Direction`).
    _Direction,

    _Normal,
    _Tangent,
    _Bitangent,
    _Barycentric,

    _Rotor,
    _Bivector,

    _Velocity,
    _Acceleration,

    _Radiance,
    _Irradiance,
    _Exitance,
    _Reflectance,
    _Absorptance,
    _Spectrum,
    // Semantic End

    // Coordinate Space States - In what basis vector lives
    _WorldSpace,     // Global scene coordinates
    _CameraSpace,    // Camera-relative
    _ClipSpace,      // After projection
    _NDCSpace,       // Normalized device coords [-1,1] or [0,1]
    _ScreenSpace,    // Pixel coordinates
    _ObjectSpace,    // Model-local
    _FrameSpace,     // Local frame space
    // Coordinate Space States End
);

pub trait ScaleSeal: Copy + Clone + 'static {
    const IS_UNIT: bool;
}
impl ScaleSeal for _Unit {
    const IS_UNIT: bool = true;
}
impl ScaleSeal for _NonUnit {
    const IS_UNIT: bool = false;
}

pub trait SpaceSeal: Copy + Clone + 'static {}
impl SpaceSeal for _WorldSpace {}
impl SpaceSeal for _CameraSpace {}
impl SpaceSeal for _ClipSpace {}
impl SpaceSeal for _NDCSpace {}
impl SpaceSeal for _ScreenSpace {}
impl SpaceSeal for _ObjectSpace {}
impl SpaceSeal for _FrameSpace {}

pub trait FreelyMutable: Copy + 'static {}

impl FreelyMutable for _Direction {}
impl FreelyMutable for _Radiance {}
impl FreelyMutable for _Irradiance {}
impl FreelyMutable for _Exitance {}
impl FreelyMutable for _Reflectance {}
impl FreelyMutable for _Absorptance {}
impl FreelyMutable for _Spectrum {}

/// Geometric Vector Marker trait
/// Enables geometric operations: dot products, cross products, angles, and distances.
/// Also prevents nonsensical operations between different types of vectors, like light and direction.
pub trait GeoSeal: Copy + 'static {}
impl GeoSeal for _Direction {}
impl GeoSeal for _Normal {}
impl GeoSeal for _Tangent {}
impl GeoSeal for _Bitangent {}
impl GeoSeal for _Rotor {}
impl GeoSeal for _Bivector {}
