use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use physarus::core::Mat3;

fn matmul_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("mat3x3");

    group.bench_function("my simd", |b| {
        b.iter_batched(
            || {
                // Setup: Generate fresh random matrices (NOT timed)
                let m1_data: [[f32; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
                let m2_data: [[f32; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));

                let m1 = Mat3::from_cols(m1_data);
                let m2 = Mat3::from_cols(m2_data);
                (m1, m2)
            },
            |(m1, m2)| {
                // Routine: This is the ONLY timed part
                black_box(&m1) * black_box(&m2)
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("glam", |b| {
        b.iter_batched(
            || {
                let m1_data: [[f32; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
                let m2_data: [[f32; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));

                let m1 = glam::Mat3A::from_cols_array_2d(&m1_data);
                let m2 = glam::Mat3A::from_cols_array_2d(&m2_data);
                (m1, m2)
            },
            |(m1, m2)| black_box(m1) * black_box(m2),
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra", |b| {
        b.iter_batched(
            || {
                let m1_data: [[f32; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
                let m2_data: [[f32; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));

                let m1 = nalgebra::Matrix3::from_fn(|i, j| m1_data[j][i]);
                let m2 = nalgebra::Matrix3::from_fn(|i, j| m2_data[j][i]);
                (m1, m2)
            },
            |(m1, m2)| black_box(&m1) * black_box(&m2),
            BatchSize::SmallInput,
        );
    });
    //
    // group.bench_function("faer", |b| {
    //     b.iter_batched(
    //         || {
    //             let m1_data: [[f32; 3]; 3] =
    //                 std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
    //             let m2_data: [[f32; 3]; 3] =
    //                 std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
    //
    //             let m1 = faer::Mat::from_fn(3, 3, |i, j| m1_data[j][i]);
    //             let m2 = faer::Mat::from_fn(3, 3, |i, j| m2_data[j][i]);
    //             (m1, m2)
    //         },
    //         |(m1, m2)| black_box(&m1) * black_box(&m2),
    //         BatchSize::SmallInput,
    //     );
    // });

    group.finish();
}

criterion_group!(benches, matmul_comparison);
criterion_main!(benches);
