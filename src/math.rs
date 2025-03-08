use rand::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

/// Returns a * b mod m
pub fn modmul(a: u64, b: u64, m: u64) -> u64 {
    let mul = (a % m) as u128 * (b % m) as u128;
    (mul % (m as u128)) as u64
}

/// Returns a + b mod m
pub fn modadd(a: u64, b: u64, m: u64) -> u64 {
    (a + b) % m
}

/// Returns a - b mod m
pub fn modsub(a: u64, b: u64, m: u64) -> u64 {
    (a + m - (b % m)) % m
}

/// Returns a^b mod m
pub fn modexp(x: u64, y: u64, p: u64) -> u64 {
    let (mut i, mut r, mut z) = (y, 1, x);
    while i != 0 {
        if i & 1 == 1 {
            r = modmul(r, z, p);
        }
        z = modmul(z, z, p);
        i >>= 1;
    }
    r
}

/// Returns a^-1 mod m
pub fn modinv(a: u64, m: u64) -> u64 {
    modexp(a, m - 2, m)
}

/// Returns the low and high bits of a * b
pub fn mul64(a: u64, b: u64) -> (u64, u64) {
    let c = a as u128 * b as u128;
    (c as u64, (c >> 64) as u64)
}

/// Returns a^-1 mod 2^64
pub fn inv(a: u64) -> u64 {
    let (mut r, mut m) = (1u64, 2u64);
    while m != 0 {
        r |= a.wrapping_mul(r) & m;
        m <<= 1;
    }
    r
}

/// Returns a primitive m-th root of unity mod p
pub fn find_primitive_root(p: u64, m: u64) -> u64 {
    let mut r;
    loop {
        r = rand::rng().random_range(2..p - 1);
        if modexp(r, (p - 1) >> 1, p) == 1 {
            break;
        }
    }
    modexp(r, (p - 1) * m, p)
}

/// Miller-Rabin prime test
fn prime_test(p: u64) -> bool {
    if p < 2 || (p != 2 && p % 2 == 0) {
        return false;
    }

    let s = p - 1;
    let x = s.leading_zeros();
    for _ in 0..256 {
        let r = rand::rng().random_range(1..p);
        let mut t = modexp(r, s, p);
        let mut j = 0;
        while j < x && t != 1 && t != p - 1 {
            t = modmul(t, t, p);
            j += 1;
        }
        if t != p - 1 && j > 1 {
            return false;
        }
    }

    true
}

/// Fills a vector with l-bit primes of the form
/// p = k * m + 1
pub fn gen_primes(l: u8, m: u64, p: &mut [u64]) {
    static K: AtomicU64 = AtomicU64::new(0);
    let lbits = 1u64 << l;
    for i in p.iter_mut() {
        loop {
            *i = lbits + K.fetch_add(1, Ordering::Relaxed) * m + 1;
            if prime_test(*i) {
                break;
            }
        }
    }
}

/// Compute the forward NTT
pub fn ntt(roots: &[u64], x: &mut [u64], n: usize, q: u64, qinv: u64) {
    let (mut m, mut t) = (1usize, n / 2);
    while m < n {
        let mut k = 0;
        for i in 0..m {
            let s = roots[m + i];
            for j in k..k + t {
                let (lo, mut hi) = mul64(x[j + t], s);
                let (_, carry) = mul64(lo.wrapping_mul(qinv), q);
                hi = if hi < carry {
                    hi + q - carry
                } else {
                    hi - carry
                };
                x[j + t] = if x[j] < hi {
                    x[j].wrapping_add(q) - hi
                } else {
                    x[j] - hi
                };
                x[j] = x[j].wrapping_add(hi);
                if x[j] > q {
                    x[j] -= q;
                }
            }
            k += t << 1;
        }
        m <<= 1;
        t >>= 1;
    }
}

/// Compute the inverse NTT
pub fn intt(iroots: &[u64], x: &mut [u64], n: usize, q: u64, qinv: u64, ninv: u64) {
    let (mut m, mut t) = (n / 2, 1usize);
    while m > 0 {
        let mut k = 0;
        for i in 0..m {
            let s = iroots[m + i];
            for j in k..k + t {
                let carry = if x[j] < x[j + t] {
                    x[j] + q - x[j + t]
                } else {
                    x[j] - x[j + t]
                };
                x[j] += x[j + t];
                if x[j] > q {
                    x[j] -= q;
                }
                let (lo, hi) = mul64(carry, s);
                let (_, carry) = mul64(lo.wrapping_mul(qinv), q);
                x[j + t] = if hi < carry {
                    hi + q - carry
                } else {
                    hi - carry
                }
            }
            k += t << 1;
        }
        m >>= 1;
        t <<= 1;
    }

    for i in x.iter_mut().take(n) {
        let (lo, hi) = mul64(*i, ninv);
        let (_, carry) = mul64(lo.wrapping_mul(qinv), q);
        *i = if hi < carry {
            hi + q - carry
        } else {
            hi - carry
        };
    }
}
