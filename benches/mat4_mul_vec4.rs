use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use physarus::core::{Direction4, Mat4};

fn mat4_transform_vec4(c: &mut Criterion) {
    let mut group = c.benchmark_group("mat4_transform_vec4");

    group.bench_function("my simd", |b| {
        b.iter_batched(
            || {
                // Generate random matrix and vector
                let m_data: [[f32; 4]; 4] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
                let v_data: [f32; 4] = std::array::from_fn(|_| rand::random());

                let m = Mat4::from_cols(m_data);
                let v = Direction4::from_array(v_data);
                (m, v)
            },
            |(m, v)| {
                // Matrix * Vector multiplication
                black_box(&m) * black_box(&v)
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("glam", |b| {
        b.iter_batched(
            || {
                let m_data: [[f32; 4]; 4] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
                let v_data: [f32; 4] = std::array::from_fn(|_| rand::random());

                let m = glam::Mat4::from_cols_array_2d(&m_data);
                let v = glam::Vec4::from_array(v_data);
                (m, v)
            },
            |(m, v)| black_box(m) * black_box(v),
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra", |b| {
        b.iter_batched(
            || {
                let m_data: [[f32; 4]; 4] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rand::random()));
                let v_data: [f32; 4] = std::array::from_fn(|_| rand::random());

                let m = nalgebra::Matrix4::from_fn(|i, j| m_data[j][i]);
                let v = nalgebra::Vector4::from_fn(|i, _| v_data[i]);
                (m, v)
            },
            |(m, v)| black_box(&m) * black_box(&v),
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, mat4_transform_vec4);
criterion_main!(benches);
