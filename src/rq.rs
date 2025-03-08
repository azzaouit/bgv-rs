use crate::math::*;
use rand::prelude::*;
use rand_distr::{Distribution, Normal};
use rug::Integer;
use std::thread;

/// Defines a polynomial ring
#[derive(Debug)]
pub struct Rq {
    /// Polynomial degree
    n: usize,
    /// Max levels supported
    l: usize,
    /// Ciphertext modulus Q
    q: Integer,
    /// Q/2
    q_half: Integer,
    /// CRT residues q_i
    qi: Vec<u64>,
    /// Q/q_i
    q_over_qi: Vec<Integer>,
    /// [Q/q_i]^-1 mod q_i
    q_over_qi_inv: Vec<u64>,
    /// [q_i]^-1
    qi_inv: Vec<u64>,
    /// [n]_{q_i}^-1
    ninv: Vec<u64>,
    /// Primitive roots of unity
    roots: Vec<u64>,
    /// Inverse primitive roots of unity
    iroots: Vec<u64>,
    /// Sampling mean
    mu: f64,
    /// Sampling standard deviation
    sigma: f64,
}

impl Rq {
    /// Initialize a polynomial ring
    pub fn new(n: usize, l: usize, nbits: u8, mu: f64, sigma: f64) -> Self {
        let nn = 2 * n as u64;
        let mut qi = vec![0u64; l];
        let mut q_over_qi: Vec<Integer> = Vec::with_capacity(l);
        let mut q_over_qi_inv: Vec<u64> = Vec::with_capacity(l);
        let mut qi_inv: Vec<u64> = Vec::with_capacity(l);
        let mut ninv: Vec<u64> = Vec::with_capacity(l);
        let mut roots = vec![0u64; l * n];
        let mut iroots = vec![0u64; l * n];

        const M32: u64 = 1u64 << 32;

        gen_primes(nbits, nn, &mut qi);

        for i in 0..l {
            let root = find_primitive_root(qi[i], nn);
            let iroot = modinv(root, qi[i]);
            assert_eq!(modexp(root, nn, qi[i]), 1);
            assert_eq!(modexp(iroot, nn, qi[i]), 1);

            qi_inv.push(inv(qi[i]));
            ninv.push(modinv(n as u64, qi[i]));
            ninv[i] = modmul(ninv[i], M32, qi[i]);
            ninv[i] = modmul(ninv[i], M32, qi[i]);

            let (mut power, mut ipower) = (1u64, 1u64);
            for j in 0..n {
                let index = i * n + ((j as u32).reverse_bits() as usize >> (32 - n.ilog2()));
                roots[index] = modmul(power, M32, qi[i]);
                roots[index] = modmul(roots[index], M32, qi[i]);
                iroots[index] = modmul(ipower, M32, qi[i]);
                iroots[index] = modmul(iroots[index], M32, qi[i]);
                power = modmul(power, root, qi[i]);
                ipower = modmul(ipower, iroot, qi[i]);
            }
        }

        let mut q = Integer::from(qi[0]);
        for i in qi.iter().take(l).skip(1) {
            q *= *i;
        }
        let q_half = q.clone() / 2;

        for i in 0..l {
            q_over_qi.push(q.clone() / qi[i]);
            q_over_qi_inv.push((q_over_qi[i].clone() % Integer::from(qi[i])).to_u64_wrapping());
            q_over_qi_inv[i] = modinv(q_over_qi_inv[i], qi[i]);
        }

        Self {
            n,
            l,
            q,
            q_half,
            qi,
            q_over_qi,
            q_over_qi_inv,
            qi_inv,
            ninv,
            roots,
            iroots,
            mu,
            sigma,
        }
    }

    /// Zero polynomial
    pub fn zero(&self) -> Vec<u64> {
        vec![0u64; self.l * self.n]
    }

    /// Encode a polynomial into CRT form
    pub fn encode(&self, x: &[u64]) -> Vec<u64> {
        let mut out = vec![0u64; self.l * self.n];
        for i in 0..self.l {
            for j in 0..self.n {
                out[i * self.n + j] = x[j] % self.qi[i];
            }
        }
        self.ntt(&mut out);
        out
    }

    /// Decode a polynomial from CRT form
    pub fn decode(&self, v: &[u64], m: u32) -> Vec<u64> {
        let mut x = v.to_owned();
        self.intt(&mut x);
        let mut out = vec![0u64; self.n];
        for i in 0..self.n {
            let mut u = Integer::new();
            for j in 0..self.l {
                let mut v = self.q_over_qi[j].clone() * self.q_over_qi_inv[j];
                v *= x[j * self.n + i];
                u += v;
            }
            u %= &self.q;
            if u > self.q_half {
                u -= &self.q;
            }
            out[i] = u.mod_u(m) as u64;
        }
        out
    }

    /// Compute the forward NTT
    pub fn ntt(&self, x: &mut [u64]) {
        for i in 0..self.l {
            thread::scope(|s| {
                s.spawn(|| {
                    let offset = i * self.n;
                    ntt(
                        &self.roots[offset..offset + self.n],
                        &mut x[offset..offset + self.n],
                        self.n,
                        self.qi[i],
                        self.qi_inv[i],
                    )
                });
            });
        }
    }

    /// Compute the inverse NTT
    pub fn intt(&self, x: &mut [u64]) {
        for i in 0..self.l {
            thread::scope(|s| {
                s.spawn(|| {
                    let offset = i * self.n;
                    intt(
                        &self.iroots[offset..offset + self.n],
                        &mut x[offset..offset + self.n],
                        self.n,
                        self.qi[i],
                        self.qi_inv[i],
                        self.ninv[i],
                    )
                });
            });
        }
    }

    /// Add a polynomial
    pub fn add_eq(&self, a: &mut [u64], b: &[u64]) {
        for i in 0..self.l {
            for j in 0..self.n {
                let index = i * self.n + j;
                a[index] = modadd(a[index], b[index], self.qi[i]);
            }
        }
    }

    /// Subtract a polynomial
    pub fn sub_eq(&self, a: &mut [u64], b: &[u64]) {
        for i in 0..self.l {
            for j in 0..self.n {
                let index = i * self.n + j;
                a[index] = modsub(b[index], a[index], self.qi[i]);
            }
        }
    }

    /// Multiply two polynomials
    pub fn mul(&self, c: &mut [u64], a: &[u64], b: &[u64]) {
        for i in 0..self.l {
            for j in 0..self.n {
                let index = i * self.n + j;
                c[index] = modmul(a[index], b[index], self.qi[i]);
            }
        }
    }

    /// Subtract a polynomial
    pub fn mul_eq(&self, a: &mut [u64], b: &[u64]) {
        for i in 0..self.l {
            for j in 0..self.n {
                let index = i * self.n + j;
                a[index] = modmul(a[index], b[index], self.qi[i]);
            }
        }
    }

    /// Multiply a polynomial by a constant
    pub fn cmul(&self, a: &mut [u64], b: u64) {
        for i in 0..self.l {
            for j in 0..self.n {
                let index = i * self.n + j;
                a[index] = modmul(a[index], b, self.qi[i]);
            }
        }
    }

    /// Sample a polynomial with uniformly random coefficients
    pub fn sample_uniform(&self) -> Vec<u64> {
        let mut x = self.zero();
        let mut rng = rand::rng();
        for i in 0..self.n {
            let s: u64 = rng.random();
            for j in 0..self.l {
                let index = j * self.n + i;
                x[index] = s % self.qi[j];
            }
        }
        x
    }

    /// Sample a polynomial with random coefficients in the range {-1, 0, 1}
    pub fn sample_ternary(&self) -> Vec<u64> {
        let mut x = self.zero();
        let mut rng = rand::rng();
        for i in 0..self.n {
            let s: i64 = rng.random_range(-1..2);
            for j in 0..self.l {
                let index = j * self.n + i;
                x[index] = (s % self.qi[j] as i64) as u64;
            }
        }
        x
    }

    /// Sample a polynomial with coefficients from a gaussian distribution
    pub fn sample_gaussian(&self) -> Vec<u64> {
        let mut x = self.zero();
        let normal = Normal::new(self.mu, self.sigma).unwrap();
        for i in 0..self.n {
            let s = normal.sample(&mut rand::rng()).floor() as i64;
            for j in 0..self.l {
                let index = j * self.n + i;
                x[index] = (s % self.qi[j] as i64) as u64;
            }
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: u32 = 65537;
    const N: usize = 1024;
    const L: usize = 15;
    const NBITS: u8 = 50;
    const MU: f64 = 0.0;
    const SIGMA: f64 = 3.19;

    #[test]
    fn test_encode() {
        let r: Rq = Rq::new(N, L, NBITS, MU, SIGMA);
        let v: Vec<u64> = (0..r.n)
            .map(|_| rand::rng().random_range(0..T) as u64)
            .collect();
        let encoded = r.encode(&v);
        let decoded = r.decode(&encoded, T);
        assert_eq!(decoded, v);
    }
}
