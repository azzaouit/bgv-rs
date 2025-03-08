# bgv-rs

An efficient FHE library

# Quickstart
```toml
[dependencies]
bgv = "0.6.0"
```

# Example

```rust
use rand::prelude::*;
use bgv::{math, BGV};

const T: u64 = 65537; // plaintext modulus
const N: usize = 1024; // polynomial degree
const L: usize = 15; // levels
const NBITS: u8 = 50; // small prime bit length
const MU: f64 = 0.0; // sampling mean
const SIGMA: f64 = 3.19; // sampling std. dev

let mut rng = rand::rng();
let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
let k = b.key_gen(); // generate keys

let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
let sum: Vec<u64> = (0..N).map(|i| math::modadd(m1[i], m2[i], T)).collect();

let x = b.encrypt(&m1, &k.pk);
let y = b.encrypt(&m2, &k.pk);
let z = b.add(&x, &y); // homomorphic addition

let m3 = b.decrypt(&z, &k.s);
assert_eq!(m3, sum);
```
