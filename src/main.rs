#![allow(
    dead_code,
    unused_imports,
    unused_mut,
    unused_variables,
    unused_must_use,
    unused_parens
)]
use std::{hint::black_box, time::Instant};

use physarus::prelude::*;

//
fn main() -> Result<(), PbrtError> {
    let mut rng = rand::rng();

    let radi = rand::random::<Radiance3>();

    let dir_v = rand::random::<Direction3>().normalize();
    let dir_u = rand::random::<Direction3>().normalize();
    let dir_w = rand::random::<Direction3>().normalize();

    let dir4 = rand::random::<Direction4>().normalize();

    let normal = Normal3::from_cross(&dir_w, &dir_v);
    let m4r1 = rand::random::<Mat4>();
    let m4r2 = rand::random::<Mat4>();
    let m4inv = m4r1.inverse().unwrap();
    // let trn1 = Transform::new(m4, m4inv);

    std::hint::black_box({
        let r = m4r1 * m4r2;
        dbg!(r);
    });

    Ok(())
}

//
// fn main() -> Result<(), PbrtError> {
//     let mut rng = rand::rng();
//
//     let radi = rand::random::<Radiance3>();
//
//     let dir_v = rand::random::<Direction3>().normalize();
//     let dir_u = rand::random::<Direction3>().normalize();
//     let dir_w = rand::random::<Direction3>().normalize();
//     //
//     // let sph_tir_area = spherical_triangle_area(dir_v, dir_u, dir_w);
//     // let sph_tir_area2 = spherical_triangle_area2(dir_v, dir_u, dir_w);
//     // dbg!(
//     //     sph_tir_area,
//     //     sph_tir_area,
//     //     sph_tir_area.eq_eps(&sph_tir_area2)
//     // );
//     //
//     // dbg!((dir_v, dir_u));
//     //
//     // let r = Rotor::from_vecs(&dir_v, &dir_u);
//     //
//     // dbg!(r.rotate(&dir_v));
//     //
//     let n = Normal3::from_cross(&dir_w, &dir_u);
//     let (n, b1, b2) = n.coordinate_system().ensure_map_err(
//         |(n, ..)| !n.is_zero_eps(),
//         |(n, ..)| {
//             pbrt_err!("One of the vectors is zero vector: {:?}", n)
//                 .context("what can i say, this sucks")
//         },
//     )?;
//     dbg!((b1.cross(&b2), &n));
//     dbg!(b1.cross(&b2).dot(&n));
//     //
//     // dbg!((n.dot(&b1), n.dot(&b2), b1.dot(&b2)));
//     // dbg!((n.len(), b1.dot(&b1), b2.dot(&b2)));
//
//     // let bb1 = Bounds3::new(Point3::ZERO, Point3::new(2.0e10, 3.0, 5.0));
//     // let bb2 = Bounds3::new(Point3::ZERO, Point3::new(6.0, 5.0, 5.0));
//     // let bb0 = Bounds3::from(Point3::ZERO);
//     // let p = Point3::new(1.0, 1.0, 1.0);
//     // dbg!(bb1, bb2);
//     // dbg!(bb1.overlaps(&bb2));
//     // dbg!(bb1.inside(&p));
//     // dbg!(bb0.distance(&p));
//     //
//     // dbg!(bb1);
//     // dbg!(bb1.max_dimension());
//
//     // let dir_v1 = Direction3::new(1.0, 1.0, 0.0).normalize();
//     // let dir_v2 = Direction3::new(0.999994, 1.0, 0.0).normalize();
//     // dbg!(dir_v1, dir_v2);
//     // dbg!(dir_v1.angle(&dir_v2));
//     // dbg!(dir_v1.dot(&dir_v2).acos());
//     // dbg!(dir_v1.angle_n(&dir_v2));
//     // dbg!(dir_v1.cross(&dir_v2).len().asin());
//     // dbg!((2.0 * (dir_v2 - dir_v1).len().atan2((dir_v1 + dir_v2).len())));
//
//     // println!("NEW NEW NEW");
//     //
//     // let dir_v = rand::random::<Direction3>().normalize();
//     // let dir_u = rand::random::<Direction3>().normalize();
//     // let dir_w = rand::random::<Direction3>().normalize();
//     //
//     // let octa = dbg!(Octahedral2::from(dir_v));
//     // let dir_v_from_octa = Unit3::from(octa);
//     //
//     // dbg!(dir_v, octa, dir_v_from_octa);
//     // dbg!(dir_v.eq_eps(&dir_v_from_octa));
//     //
//     // let p = Point2::new(0.5, 0.5);
//     // dbg!(p);
//     // let sph_v = equal_area_square_to_sphere(p);
//     // let p = equal_area_sphere_to_square(&sph_v);
//     //
//     // dbg!(sph_v, p);
//     //
//     // let x_axis = unsafe { Unit3::neu([1.0, 0.0, 0.0]) };
//     // let y_axis = unsafe { Unit3::neu([0.0, 1.0, 0.0]) };
//     // let z_axis = unsafe { Unit3::neu([0.0, 0.0, 1.0]) };
//     // let rotor = Rotor::from_axis_angle(&(x_axis + y_axis).normalize(), Rad::new(PI / 4.0));
//     // let rotor2 = Rotor::from_vecs_angle(&x_axis, &y_axis, Rad::new(PI / 4.0));
//     // let v = Direction3::new(1.0, 0.0, 0.0);
//     // let rotor3 = Rotor::from_vecs(&x_axis, &y_axis);
//     // dbg!("AFTER HEHE");
//     // dbg!(v);
//     // dbg!(rotor.rotate(&v));
//     // dbg!(rotor2.rotate(&v));
//     // dbg!(rotor3.rotate(&v));
//     // dbg!(rotor3);
//     //
//     // let v2d1 = rand::random::<Direction2>().normalize();
//     // let v2d2 = rand::random::<Direction2>().normalize();
//     //
//     // dbg!(v2d1.angle_n(&v2d2));
//     // dbg!(v2d1.angle(&v2d2));
//     //
//     // let dir_v = rand::random::<Direction3>();
//     // let dir_u = rand::random::<Direction3>();
//     // let dir_w = rand::random::<Direction3>();
//     // dbg!(dir_v.angle_n(&dir_w));
//     // dbg!(dir_v.angle(&dir_w));
//     //
//     let v1 = Direction3::new(33962.035, 41563.4, 7706.415);
//     let v2 = Direction3::new(-24871.969, -30438.8, -5643.727);
//
//     v1.normalize().len_sqr();
//
//     let v1crossv2 = v1.cross(&v2);
//
//     dbg!(v1crossv2);
//     dbg!(v1crossv2.is_perpendicular(&v1));
//     dbg!(v1crossv2.is_perpendicular(&v2));
//     dbg!(v1crossv2.dot(&v1));
//     dbg!(v1crossv2.dot(&v2));
//
//     dbg!(v1.dot(&v2));
//
//     let t = Instant::now();
//     // for i in 0..10000 {
//     // let m1 = rand::random::<Mat3>();
//     //
//     // let m1_inv = m1.inverse().unwrap();
//     // dbg!(m1 * m1_inv);
//     // }
//
//     let normal = Normal3::from_cross(&dir_w, &dir_v);
//     let dir = rand::random::<Direction3>();
//
//     let m4 = rand::random::<Mat3>().to_hom_4x4() * 100.0;
//     let m4inv = m4.inverse().unwrap();
//
//     let trn1 = Transform::new(m4, m4inv);
//     let m4 = rand::random::<Mat3>().to_hom_4x4() * 100.0;
//     let m4inv = m4.inverse().unwrap();
//     let trn2 = Transform::new(m4, m4inv);
//     let m4 = rand::random::<Mat3>().to_hom_4x4() * 100.0;
//     let m4inv = m4.inverse().unwrap();
//     let trn3 = Transform::new(m4, m4inv);
//
//     std::hint::black_box({
//         let tn = normal.apply_transform(&trn1);
//         dbg!(tn);
//     });
//
//     let xxzy = swizzle!(dir => x, x, z, y);
//
//     dbg!(dir);
//
//     let angle = PI * rand::random::<Float>();
//     let rotor = Rotor::from_axis_angle(&dir.normalize(), Rad::new(angle));
//     dbg!(rotor.rotate(&Direction3::Z_AXIS));
//
//     let tr = Transform::from_rotor(&rotor);
//     dbg!(tr.transform(&Direction3::Z_AXIS));
//
//     let gmat4 = glam::Mat4::from_axis_angle(
//         glam::Vec3::from_array([dir.x(), dir.y(), dir.z()]).normalize(),
//         angle,
//     );
//
//     dbg!(gmat4.transform_vector3(glam::Vec3::Z));
//
//     Ok(())
// }
