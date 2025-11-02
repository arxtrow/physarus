use crate::prelude::*;

pub fn spherical_triangle_area(a: Unit3, b: Unit3, c: Unit3) -> Float {
    let det = a.dot(&b.cross(&c));
    let den = 1.0 + a.dot(&b) + a.dot(&c) + b.dot(&c);

    (2.0 * det.atan2(den)).abs()
}

pub fn spherical_triangle_area2(a: Unit3, b: Unit3, c: Unit3) -> Float {
    let axb = a.cross(&b);
    let bxc = b.cross(&c);
    let cxa = c.cross(&a);

    if [axb, bxc, cxa].iter().any(|v| v.len_sqr().is_zero_eps()) {
        return 0.0;
    }

    let bxa = -axb;
    let cxb = -bxc;
    let axc = -cxa;

    let Rad(alpha) = cxa.angle(&bxa);
    let Rad(betta) = axb.angle(&cxb);
    let Rad(gamma) = bxc.angle(&axc);

    (alpha + betta + gamma - PI).abs()
}

pub fn spherical_quad_area(a: Unit3, b: Unit3, c: Unit3, d: Unit3) -> Float {
    let axb = a.cross(&b);
    let bxc = b.cross(&c);
    let cxd = c.cross(&d);
    let dxa = d.cross(&a);

    if [axb, bxc, cxd, dxa]
        .iter()
        .any(|v| v.len_sqr().is_zero_eps())
    {
        return 0.0;
    }

    let bxa = -axb;
    let cxb = -bxc;
    let dxc = -cxd;
    let axd = -dxa;

    let Rad(alpha) = dxa.angle(&bxa);
    let Rad(betta) = axb.angle(&cxb);
    let Rad(gamma) = bxc.angle(&dxc);
    let Rad(delta) = cxd.angle(&axd);

    (alpha + betta + gamma + delta - 2.0 * PI).abs()
}

pub fn spherical_direction(sin_theta: Float, cos_theta: Float, phi: Float) -> Unit3 {
    // SAFETY<Unit> : cos^2(𝜑)sin^2(𝜃) + sin^2(𝜑)sin^2(𝜃) + cos^2(𝜃)
    //                = sin^2(𝜃)(cos^2(𝜑) + sin^2(𝜑)) + cos^2(𝜃) = sin^2(𝜃) * 1 + cos^2(𝜃) = 1
    unsafe {
        Vec3::raw([
            phi.cos() * sin_theta.clamp(-1.0, 1.0),
            phi.sin() * sin_theta.clamp(-1.0, 1.0),
            cos_theta.clamp(-1.0, 1.0),
        ])
    }
}

pub fn spherical_theta(v: Unit3) -> Rad {
    Rad::new(v.z().acos())
}

pub fn spherical_phi(v: Unit3) -> Rad {
    let p = v.y().atan2(v.x());

    Rad::new(if p < 0.0 { p + 2.0 * PI } else { p })
}

/// only if n aligned with z-axis
pub fn cos_theta(v: Unit3) -> Float {
    v.z()
}

pub fn cos_sqr_theta(v: Unit3) -> Float {
    v.z() * v.z()
}

pub fn sin_theta(v: Unit3) -> Float {
    sin_sqr_theta(v).sqrt()
}

pub fn sin_sqr_theta(v: Unit3) -> Float {
    (1.0 - cos_sqr_theta(v)).max(0.0)
}

pub fn tan_theta(v: Unit3) -> Float {
    sin_theta(v) / cos_theta(v)
}

pub fn tan_theta_sqr(v: Unit3) -> Float {
    sin_sqr_theta(v) / cos_sqr_theta(v)
}

pub fn cos_phi(v: Unit3) -> Float {
    let radius = sin_theta(v);

    if radius.is_zero_eps() {
        return 1.0;
    }

    (v.x() / radius).clamp(-1.0, 1.0)
}

pub fn sin_phi(v: Unit3) -> Float {
    let radius = sin_theta(v);

    if radius.is_zero_eps() {
        return 0.0;
    }

    (v.y() / radius).clamp(-1.0, 1.0)
}

pub fn cos_delta_phi(va: Unit3, vb: Unit3) -> Float {
    let plane_dot = va.x() * vb.x() + va.y() * vb.y();

    let plane_va_len_sqr = va.x() * va.x() + va.y() * va.y();
    let plane_vb_len_sqr = vb.x() * vb.x() + vb.y() * vb.y();

    if plane_va_len_sqr.is_zero_eps() || plane_vb_len_sqr.is_zero_eps() {
        return 1.0;
    }

    (plane_dot / (plane_va_len_sqr * plane_vb_len_sqr).sqrt()).clamp(-1.0, 1.0)
}

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Octahedral2 {
    x: u16,
    y: u16,
}

impl From<Unit3> for Octahedral2 {
    fn from(v: Unit3) -> Self {
        let v = v / (v.x().abs() + v.y().abs() + v.z().abs());

        if v.z() >= 0.0 {
            let x = encode_float_to_u16(v.x());
            let y = encode_float_to_u16(v.y());

            Self { x, y }
        } else {
            // Reflect across x + y = 1 and restore quadrant by signum
            let x = (1.0 - v.y().abs()) * v.x().signum();
            let y = (1.0 - v.x().abs()) * v.y().signum();

            let x = encode_float_to_u16(x);
            let y = encode_float_to_u16(y);

            Self { x, y }
        }
    }
}

impl From<Octahedral2> for Unit3 {
    fn from(v: Octahedral2) -> Self {
        let mut x = decode_u16_to_float(v.x);
        let mut y = decode_u16_to_float(v.y);
        let z = 1.0 - (x.abs() + y.abs());

        if z < 0.0 {
            let x_reflected = x;
            let y_reflected = y;

            x = (1.0 - y_reflected.abs()) * x_reflected.signum();
            y = (1.0 - x_reflected.abs()) * y_reflected.signum();
        }

        Direction3::new(x, y, z).normalize()
    }
}

/// Performs the encoding from a value in `[-1.0, 1.0]` to the `[0, u16::MAX]`.
///
/// Note: Due to quantization, `decode(encode(x))` may differ from `x` by up to
/// `1.0 / 65535.0 ≈ 0.0000153` (half the quantization step).
fn encode_float_to_u16(x: Float) -> u16 {
    debug_assert!((-1.0..1.0).contains(&x));

    let x_0_to_1 = ((x + 1.0) / 2.0).clamp(0.0, 1.0);

    (x_0_to_1 * (u16::MAX as Float)).round() as u16
}

/// Performs the decoding from a value in `[0, u16::MAX]` to the [-1.0, 1.0] float.
fn decode_u16_to_float(x: u16) -> Float {
    let x_0_to_1 = x as Float / u16::MAX as Float;

    (x_0_to_1 * 2.0) - 1.0
}

pub fn equal_area_square_to_sphere(p: Point2) -> Unit3 {
    debug_assert!((0.0..=1.0).contains(&p.x()) && (0.0..=1.0).contains(&p.y()));

    // Transform p from [0, 1]² to [-1, 1]²
    let u = p.x() * 2.0 - 1.0;
    let v = p.y() * 2.0 - 1.0;

    // Work in first quadrant, restore signs at end
    let u_abs = u.abs();
    let v_abs = v.abs();

    // Decompose point parallel/perpendicular to diagonal in first quadrant
    let a = u_abs + v_abs; // == |diag.proj(p)|
    let b = v_abs - u_abs; // == |diag.reject(p)|

    // distance from point to l: x + y = 1.0 computed in first quadrant using |u| and |v|
    // unnormalized; true Euclidean would be (1.0 - |u| + |v|) / sqrt(2.0))
    // but constant factor doesn't matter for the concentric mapping(??)
    let signed_dist = 1.0 - a;
    let dist = signed_dist.abs();

    // Northern hemisphere: |u| + |v| <= 1 => r_disk = u + v => r_disk in [0.0, 1.0]
    // Southern hemisphere: |u| + |v| > 1 => r_disk = 2 - (u + v) => r_disk in [0.0, 1.0]
    let r_disk = 1.0 - dist;

    let phi = if r_disk.is_zero_eps() {
        FRAC_PI_4 // Undefined at poles; pi/4 ensures continuity
    } else {
        // as phi goes from 0 to pi/4 to pi/2, b/a goes from -1.0 to 0.0 to 1.0
        // => phi = k(b/a + 1.0)
        // at (0, 1): b/a = 1 => phi = 2k and phi = pi/2 at (0, 1) => k = pi/4
        // phi = pi/4(b/a + 1)
        let a = r_disk; // to include southern hemisphere
        (b / a + 1.0) * (FRAC_PI_4)
    };

    // if signed_dist >= 0 => |u| + |v| <= 1.0 => northern hemisphere
    // if signed_dist < 0 => |u| + |v| > 1.0 => southern hemisphere
    let z = (1.0 - r_disk * r_disk).copysign(signed_dist);

    // Compute cos and sin and restore original quadrant
    let cos_phi = phi.cos().copysign(u);
    let sin_phi = phi.sin().copysign(v);

    // SAFETY<Unit>: x² + y² + z² = [cos²φ + sin²φ]·r²(2-r²) + (1-r²)²
    //                            = r²(2-r²) + 1 - 2r² + r⁴
    //                            = 2r² - r⁴ + 1 - 2r² + r⁴ = 1
    unsafe {
        Unit3::raw([
            cos_phi * r_disk * (2.0 - r_disk * r_disk).sqrt(),
            sin_phi * r_disk * (2.0 - r_disk * r_disk).sqrt(),
            z,
        ])
    }
}

pub fn equal_area_sphere_to_square(v: &Unit3) -> Point2 {
    // Work in first quadrant, restore signs at end
    let x = v.x().abs();
    let y = v.y().abs();
    let z = v.z().abs();

    // Compute disk radius from z-coordinate
    // Forward map: z = 1 - r²
    // r = sqrt(1 - z)
    let r_disk = (1.0 - z).sqrt();

    // Compute angle φ in disk coordinates
    let phi = if r_disk.is_zero_eps() {
        // At poles, angle is undefined; use arbitrary value
        FRAC_PI_4
    } else {
        y.atan2(x)
    };

    // Convert disk (r, φ) back to square (u, v) using inverse concentric map
    // Forward map: r = u+v and φ = (π/4)[(v-u)/r + 1]
    // v-u = r(4φ/π - 1) and u+v = r
    // v = r(4φ/π - 1) + r - v => 2v = r(4φ/π)
    // v = 2r*φ/π, u = r - v (for northern hemisphere)
    let v_abs = 2.0 * phi * r_disk / PI;
    let u_abs = r_disk - v_abs;

    // Handle hemisphere
    let (u_abs, v_abs) = if v.z() < 0.0 {
        // Southern hemisphere: points are in corners, need reflection
        let u_new = 1.0 - v_abs;
        let v_new = 1.0 - u_abs;

        (u_new, v_new)
    } else {
        // Northern hemisphere
        (u_abs, v_abs)
    };

    // Restore original quadrant signs
    let u = u_abs.copysign(v.x());
    let v = v_abs.copysign(v.y());

    // Transform from [-1, 1]² to [0, 1]²
    Point2::new((u + 1.0) * 0.5, (v + 1.0) * 0.5)
}

/// Wraps points outside `[0,1]²` square back inside using mirror symmetry.
/// Preserves equal-area property for tiled sphere sampling.
pub fn wrap_equal_area_square(p: &Point2) -> Point2 {
    let mut x = p.x();
    let mut y = p.y();

    // Wrap horizontal (u-axis): mirror across edges and center
    if x < 0.0 {
        x = -x; // Reflect across u = 0
        y = 1.0 - y; // Reflect across v = 0.5
    } else if x > 1.0 {
        x = 2.0 - x; // Reflect across u = 1
        y = 1.0 - y; // Reflect across v = 0.5
    }

    // Wrap vertical (v-axis): mirror across edges and center
    if y < 0.0 {
        x = 1.0 - x; // Reflect across u = 0.5
        y = -y; // Reflect across v = 0
    } else if y > 1.0 {
        x = 1.0 - x; // Reflect across u = 0.5
        y = 2.0 - y; // Reflect across v = 1
    }

    Point2::new(x, y)
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct DirectionCone<Space: SpaceSeal> {
    d: Unit3<Space>,
    /// cosine of the cone's half-angle
    cos_theta: Float,
}

impl<Space: SpaceSeal> DirectionCone<Space> {
    pub fn new(d: Unit3<Space>, cos_theta: Float) -> Self {
        debug_assert!((-1.0..=1.0).contains(&cos_theta));

        Self { d, cos_theta }
    }

    pub fn from_half_angle(d: Unit3<Space>, half_angle: Rad) -> Self {
        Self {
            d,
            cos_theta: half_angle.cos(),
        }
    }

    pub fn containing_bounds(b: &Bounds3<Space>, from: &Point3<Space>) -> Self {
        let (p_center, radius) = b.bounding_sphere();

        if from.distance_sqr(&p_center) <= radius * radius {
            return Self::entire_sphere();
        }

        let d = (p_center - from).normalize();

        let sin_sqr_theta = radius * radius / from.distance_sqr(&p_center);
        let cos_theta = (1.0 - sin_sqr_theta).max(0.0).sqrt();

        Self { d, cos_theta }
    }

    pub fn entire_sphere() -> Self {
        Self {
            d: unsafe { Unit3::raw([0.0, 0.0, 1.0]) },
            cos_theta: -1.0,
        }
    }

    pub fn union(&self, rhs: &Self) -> Self {
        let Rad(self_theta) = self.theta();
        let Rad(rhs_theta) = rhs.theta();

        let Rad(between_theta) = self.d.angle(&rhs.d);

        // Handle the cases where one cone is inside the other
        if (between_theta + rhs_theta).min(PI) <= self_theta {
            return *self;
        }
        if (between_theta + self_theta).min(PI) <= rhs_theta {
            return *rhs;
        }

        // Neither contains the other - merge them
        let merged_theta = (self_theta + rhs_theta + between_theta) / 2.0;
        if merged_theta >= PI {
            return Self::entire_sphere();
        }

        // Compute the rotation angle from self.d to the merged direction
        let r_theta = merged_theta - self_theta;
        let d = Rotor::from_vecs_angle(&self.d, &rhs.d, Rad::new(r_theta)).rotate(&self.d);

        Self {
            d,
            cos_theta: merged_theta.cos(),
        }
    }

    pub fn is_inside(&self, v: &Unit3<Space>) -> bool {
        self.d.dot(v) >= self.cos_theta
    }

    pub fn theta(&self) -> Rad {
        Rad::new(self.cos_theta.acos())
    }
}
