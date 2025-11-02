use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use physarus::core::mat4::Mat4;

fn mat4_determinant_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("mat4_determinant");

    group.bench_function("my simd", |b| {
        b.iter_batched(
            || {
                // Setup: Generate fresh random matrix (NOT timed)
                let m_data: [[f32; 4]; 4] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));

                Mat4::from_cols(m_data)
            },
            |m| {
                // Routine: This is the ONLY timed part
                black_box(black_box(&m).det())
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("glam", |b| {
        b.iter_batched(
            || {
                let m_data: [[f32; 4]; 4] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));

                glam::Mat4::from_cols_array_2d(&m_data)
            },
            |m| black_box(black_box(m).determinant()),
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra", |b| {
        b.iter_batched(
            || {
                let m_data: [[f32; 4]; 4] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));

                nalgebra::Matrix4::from_fn(|i, j| m_data[j][i])
            },
            |m| black_box(black_box(&m).determinant()),
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, mat4_determinant_comparison);
criterion_main!(benches);
