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

fn main() -> Result<(), PbrtError> {
    let mut rng = rand::rng();

    let radi = rand::random::<Radiance3>();

    let dir_v = rand::random::<Direction3>().normalize();
    let dir_u = rand::random::<Direction3>().normalize();
    let dir_w = rand::random::<Direction3>().normalize();
    let n = Normal3::from_cross(&dir_w, &dir_u);

    Ok(())
}
