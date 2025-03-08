//! An efficient FHE library
//!
//! # Example
//!```rust
//! use rand::prelude::*;
//! use bgv::{math, BGV};
//!
//! const T: u64 = 65537;
//! const N: usize = 1024;
//! const L: usize = 15;
//! const NBITS: u8 = 50;
//! const MU: f64 = 0.0;
//! const SIGMA: f64 = 3.19;
//!
//! let mut rng = rand::rng();
//! let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
//! let k = b.key_gen();
//!
//! let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
//! let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
//! let sum: Vec<u64> = (0..N).map(|i| math::modadd(m1[i], m2[i], T)).collect();
//!
//! let x = b.encrypt(&m1, &k.pk);
//! let y = b.encrypt(&m2, &k.pk);
//! let c = b.add(&x, &y);
//!
//! let m3 = b.decrypt(&c, &k.s);
//! assert_eq!(m3, sum);
//!```

pub mod math;
pub mod rq;
use serde::{Deserialize, Serialize};

use crate::rq::Rq;

/// Implements the BGV scheme
pub struct BGV {
    /// Polynomial ring
    pub r: Rq,
    /// Plaintext modulus
    pub t: u64,
}

/// Generic Key pair
#[derive(Serialize, Deserialize)]
pub struct KeyPair(pub Vec<u64>, pub Vec<u64>);

#[derive(Serialize, Deserialize)]
pub struct Ciphertext(pub Vec<u64>, pub Vec<u64>);

/// Represents a BGV key pair
#[derive(Serialize, Deserialize)]
pub struct BGVKey {
    /// Secret key
    pub s: Vec<u64>,
    /// Public key pair
    pub pk: KeyPair,
    /// Evaluation key pair
    pub ek: KeyPair,
}

impl BGV {
    /// Initialize the BGV scheme
    pub fn new(n: usize, l: usize, nbits: u8, mu: f64, sigma: f64, t: u64) -> Self {
        let r = Rq::new(n, l, nbits, mu, sigma);
        Self { r, t }
    }

    /// Generate an evaluation key pair from a secret key
    fn eval_gen(&self, s: &[u64]) -> KeyPair {
        let mut s2 = s.to_owned();
        self.r.mul_eq(&mut s2, s);

        let mut e = self.r.sample_gaussian();
        self.r.cmul(&mut e, self.t);
        self.r.ntt(&mut e);

        let mut a = self.r.sample_uniform();
        self.r.ntt(&mut a);

        let mut b = a.clone();
        self.r.mul_eq(&mut b, s);
        self.r.sub_eq(&mut b, &e);
        self.r.add_eq(&mut b, &s2);

        KeyPair(a, b)
    }

    /// Generate a fresh key pair
    pub fn key_gen(&self) -> BGVKey {
        let mut s = self.r.sample_ternary();
        self.r.ntt(&mut s);

        let mut a = self.r.sample_uniform();
        self.r.ntt(&mut a);

        let mut e = self.r.sample_gaussian();
        self.r.cmul(&mut e, self.t);
        self.r.ntt(&mut e);

        let mut b = a.clone();
        self.r.mul_eq(&mut b, &s);
        self.r.sub_eq(&mut b, &e);

        let eval = self.eval_gen(&s);

        BGVKey {
            s,
            pk: KeyPair(a, b),
            ek: eval,
        }
    }

    /// Encrypt a polynomial
    pub fn encrypt(&self, p: &[u64], pk: &KeyPair) -> Ciphertext {
        let m = self.r.encode(p);

        let mut u = self.r.sample_ternary();
        self.r.ntt(&mut u);

        let mut e1 = self.r.sample_gaussian();
        self.r.cmul(&mut e1, self.t);
        self.r.ntt(&mut e1);

        let mut e2 = self.r.sample_gaussian();
        self.r.cmul(&mut e2, self.t);
        self.r.ntt(&mut e2);

        let mut c1 = u.clone();
        self.r.mul_eq(&mut c1, &pk.0);
        self.r.add_eq(&mut c1, &e1);

        let mut c0 = u.clone();
        self.r.mul_eq(&mut c0, &pk.1);
        self.r.add_eq(&mut c0, &e2);
        self.r.add_eq(&mut c0, &m);

        Ciphertext(c0, c1)
    }

    /// Decrypt a ciphertext
    pub fn decrypt(&self, c: &Ciphertext, s: &[u64]) -> Vec<u64> {
        let mut m = c.1.to_owned();
        self.r.mul_eq(&mut m, s);
        self.r.add_eq(&mut m, &c.0);
        self.r.decode(&m, self.t as u32)
    }

    /// Add two ciphertexts
    pub fn add(&self, a: &Ciphertext, b: &Ciphertext) -> Ciphertext {
        let mut c0 = a.0.to_owned();
        let mut c1 = a.1.to_owned();
        self.r.add_eq(&mut c0, &b.0);
        self.r.add_eq(&mut c1, &b.1);
        Ciphertext(c0, c1)
    }

    /// Multiply two ciphertexts
    pub fn mul(&self, a: &Ciphertext, b: &Ciphertext, ek: &KeyPair) -> Ciphertext {
        let mut c0 = a.0.to_owned();
        self.r.mul_eq(&mut c0, &b.0);

        let mut c2 = a.1.to_owned();
        self.r.mul_eq(&mut c2, &b.1);

        let mut c1 = a.0.to_owned();
        self.r.mul_eq(&mut c1, &b.1);

        let mut tmp = a.1.to_owned();
        self.r.mul_eq(&mut tmp, &b.0);

        self.r.add_eq(&mut c1, &tmp);

        /* Relinearize */
        let mut tmp = c2.clone();
        self.r.mul_eq(&mut tmp, &ek.1);
        self.r.mul_eq(&mut c2, &ek.0);
        self.r.add_eq(&mut c0, &tmp);
        self.r.add_eq(&mut c1, &c2);

        Ciphertext(c0, c1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math;
    use rand::prelude::*;

    const T: u64 = 65537;
    const N: usize = 1024;
    const L: usize = 15;
    const NBITS: u8 = 50;
    const MU: f64 = 0.0;
    const SIGMA: f64 = 3.19;

    #[test]
    fn test_encrypt() {
        let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
        let k = b.key_gen();
        let mut rng = rand::rng();

        let m: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();

        let c = b.encrypt(&m, &k.pk);
        let m2 = b.decrypt(&c, &k.s);

        assert_eq!(m, m2);
    }

    #[test]
    fn test_add() {
        let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
        let k = b.key_gen();
        let mut rng = rand::rng();

        let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
        let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
        let sum: Vec<u64> = (0..N).map(|i| math::modadd(m1[i], m2[i], T)).collect();

        let x = b.encrypt(&m1, &k.pk);
        let y = b.encrypt(&m2, &k.pk);
        let z = b.add(&x, &y);

        let m3 = b.decrypt(&z, &k.s);
        assert_eq!(m3, sum);
    }

    #[test]
    fn test_commutative_add() {
        let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
        let k = b.key_gen();
        let mut rng = rand::rng();

        let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
        let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();

        let x = b.encrypt(&m1, &k.pk);
        let y = b.encrypt(&m2, &k.pk);

        let z = b.add(&x, &y);
        let u = b.add(&y, &x);

        let m3 = b.decrypt(&z, &k.s);
        let m4 = b.decrypt(&u, &k.s);
        assert_eq!(m3, m4);
    }

    #[test]
    fn test_mul() {
        let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
        let k = b.key_gen();
        let mut rng = rand::rng();

        let mut m1 = vec![0u64; N];
        m1[0] = rng.random_range(0..T) as u64;

        let mut m2 = vec![0u64; N];
        m2[0] = rng.random_range(0..T) as u64;

        let x = b.encrypt(&m1, &k.pk);
        let y = b.encrypt(&m2, &k.pk);
        let z = b.mul(&x, &y, &k.ek);

        let m3 = b.decrypt(&z, &k.s);
        assert_eq!(math::modmul(m1[0], m2[0], T), m3[0]);
    }

    #[test]
    fn test_commutative_mul() {
        let b = BGV::new(N, L, NBITS, MU, SIGMA, T);
        let k = b.key_gen();
        let mut rng = rand::rng();

        let m1: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();
        let m2: Vec<u64> = (0..N).map(|_| rng.random_range(0..T) as u64).collect();

        let x = b.encrypt(&m1, &k.pk);
        let y = b.encrypt(&m2, &k.pk);

        let z = b.mul(&x, &y, &k.ek);
        let u = b.mul(&y, &x, &k.ek);

        let m3 = b.decrypt(&z, &k.s);
        let m4 = b.decrypt(&u, &k.s);
        assert_eq!(m3, m4);
    }
}
