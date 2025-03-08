use bgv::BGV;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::prelude::*;

const T: u64 = 65537;
const N: usize = 1024;
const L: usize = 15;
const NBITS: u8 = 50;
const MU: f64 = 0.0;
const SIGMA: f64 = 3.19;

fn new_benchmark(c: &mut Criterion) {
    c.bench_function("new", |b| {
        b.iter(|| {
            BGV::new(
                black_box(N),
                black_box(L),
                black_box(NBITS),
                black_box(MU),
                black_box(SIGMA),
                black_box(T),
            )
        })
    });
}

fn key_gen_benchmark(c: &mut Criterion) {
    let bgv = BGV::new(N, L, NBITS, MU, SIGMA, T);
    c.bench_function("key_gen", |b| b.iter(|| bgv.key_gen()));
}

fn encrypt_benchmark(c: &mut Criterion) {
    let mut rng = rand::rng();
    let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
    let bgv = BGV::new(N, L, NBITS, MU, SIGMA, T);
    let k = bgv.key_gen();

    c.bench_function("encrypt", |b| {
        b.iter(|| bgv.encrypt(black_box(&m1), black_box(&k.pk)))
    });
}

fn add_benchmark(c: &mut Criterion) {
    let mut rng = rand::rng();
    let bgv = BGV::new(N, L, NBITS, MU, SIGMA, T);
    let k = bgv.key_gen();
    let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
    let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
    let ct1 = bgv.encrypt(&m1, &k.pk);
    let ct2 = bgv.encrypt(&m2, &k.pk);

    c.bench_function("add", |b| {
        b.iter(|| bgv.add(black_box(&ct1), black_box(&ct2)))
    });
}

fn mul_benchmark(c: &mut Criterion) {
    let mut rng = rand::rng();
    let bgv = BGV::new(N, L, NBITS, MU, SIGMA, T);
    let k = bgv.key_gen();
    let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
    let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
    let ct1 = bgv.encrypt(&m1, &k.pk);
    let ct2 = bgv.encrypt(&m2, &k.pk);

    c.bench_function("mul", |b| {
        b.iter(|| bgv.mul(black_box(&ct1), black_box(&ct2), black_box(&k.ek)))
    });
}

criterion_group!(
    benches,
    new_benchmark,
    key_gen_benchmark,
    encrypt_benchmark,
    add_benchmark,
    mul_benchmark
);
criterion_main!(benches);
