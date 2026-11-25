//! Baillie-PSW primality test — all documented variants.
//!
//! Implements 9 BPSW variants using num-bigint for arbitrary-precision arithmetic.
//! Pure Rust, no external crypto dependencies.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

// ---------------------------------------------------------------------------
// Core mathematical primitives
// ---------------------------------------------------------------------------

/// Normalize x mod n to [0, n-1] (Rust's % can return negative).
fn mod_n(x: &BigInt, n: &BigInt) -> BigInt {
    ((x % n) + n) % n
}

/// Compute the Jacobi symbol (a/n). n must be odd and positive.
pub fn jacobi_symbol(a: &BigInt, n: &BigInt) -> i32 {
    if n <= &BigInt::zero() || n.is_even() {
        panic!("n must be odd and positive");
    }
    // Normalize a to [0, n-1] (Rust's % can return negative)
    let mut a = mod_n(a, n);
    if a.is_zero() {
        return 0;
    }
    let mut result: i32 = 1;
    let mut m = n.clone();

    while !a.is_zero() {
        while a.is_even() {
            a = a / 2;
            let r = &m % 8;
            if r == BigInt::from(3) || r == BigInt::from(5) {
                result = -result;
        }
        }
        std::mem::swap(&mut a, &mut m);
        let a_mod_4 = &a % 4;
        let m_mod_4 = &m % 4;
        if a_mod_4 == BigInt::from(3) && m_mod_4 == BigInt::from(3) {
            result = -result;
        }
        a = a % &m;
    }

    if m == BigInt::one() {
        result
    } else {
        0
    }
}

/// Check if n is a perfect square.
pub fn is_square(n: &BigInt) -> bool {
    if n < &BigInt::zero() {
        return false;
    }
    // Integer square root via Newton's method
    let mut x = n.clone();
    loop {
        let y = (x.clone() + n / &x) / 2;
        if y >= x {
            break;
        }
        x = y;
    }
    &x * &x == *n
}

/// Strong probable prime test base 2 (Miller-Rabin).
pub fn is_strong_probable_prime_base2(n: &BigInt) -> bool {
    if n < &BigInt::from(2) {
        return false;
    }
    if n == &BigInt::from(2) {
        return true;
    }
    if n.is_even() {
        return false;
    }

    // Write n - 1 = d * 2^s with d odd
    let mut d: BigInt = n - 1;
    let mut s = 0u32;
    while d.is_even() {
        d = d / 2;
        s += 1;
    }

    let x = modpow(&BigInt::from(2), &d, n);
    if x == BigInt::one() || x == n - 1 {
        return true;
    }

    let mut x = x;
    for _ in 0..s - 1 {
        x = modpow(&x, &BigInt::from(2), n); // x = x^2 mod n
        if x == n - 1 {
            return true;
        }
    }
    false
}

/// Modular exponentiation: base^exp mod modulus.
fn modpow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    base.modpow(exp, modulus)
}

/// Compute (U_k mod n, V_k mod n, Q^k mod n) using binary ladder.
/// Lucas sequences: U_0=0, U_1=1; V_0=2, V_1=P.
pub fn lucas_sequence(n: &BigInt, p: &BigInt, q: &BigInt, k: &BigInt) -> (BigInt, BigInt, BigInt) {
    if k.is_zero() {
        return (
            BigInt::zero() % n,
            BigInt::from(2) % n,
            BigInt::one() % n,
        );
    }

    let inv2: BigInt = (n + 1) / 2; // modular inverse of 2 mod n (n is odd)
    let d: BigInt = p * p - BigInt::from(4) * q; // discriminant
    let d_mod = &d % n;

    // Binary representation of k, MSB first
    let mut bits: Vec<u32> = Vec::new();
    let mut kk = k.clone();
    while !kk.is_zero() {
        bits.push(kk.to_u32().unwrap() & 1);
        kk = kk / 2;
    }
    bits.reverse();

    let mut u = BigInt::zero();   // U_0
    let mut v = BigInt::from(2);  // V_0
    let mut qk = BigInt::one();   // Q^0

    for bit in bits {
        // Double: m -> 2m
        u = mod_n(&(u.clone() * v.clone()), n);
        v = mod_n(&(v.clone() * v.clone() - 2 * qk.clone()), n);
        qk = mod_n(&(qk.clone() * qk.clone()), n);

        if bit == 1u32 {
            // Advance: 2m -> 2m+1
            let u_new = mod_n(&((p.clone() * u.clone() + v.clone()) * inv2.clone()), n);
            let v_new = mod_n(&((d_mod.clone() * u.clone() + p.clone() * v.clone()) * inv2.clone()), n);
            let qk_new = mod_n(&(qk.clone() * q), n);
            u = u_new;
            v = v_new;
            qk = qk_new;
        }
    }

    (u, v, qk)
}

/// Selfridge Method A: scan D in 5, -7, 9, -11, ... until Jacobi(D,n) = -1.
/// Returns (D, P=1, Q=(1-D)/4). Returns (D, 0, 0) if factor found.
pub fn selfridge_method_a(n: &BigInt) -> (BigInt, BigInt, BigInt) {
    let mut d = BigInt::from(5);
    for _ in 0..100000 {
        let j = jacobi_symbol(&d, n);
        if j == -1 {
            let p = BigInt::one();
            let q = (BigInt::one() - &d) / 4;
            return (d, p, q);
        } else if j == 0 {
            let g = d.abs().gcd(n);
            if g > BigInt::one() && g < *n {
                return (d, BigInt::zero(), BigInt::zero());
            }
            if d.abs() < *n {
                return (d, BigInt::zero(), BigInt::zero());
            }
        }
        // Next D: 5, -7, 9, -11, ...
        if d > BigInt::zero() {
            d = -(&d + BigInt::from(2));
        } else {
            d = -d + 2;
        }
    }
    panic!("Selfridge Method A: no suitable D found within limit");
}

/// Selfridge Method A*: same as A, but if Q == -1 (D == 5), set P = Q = 5.
pub fn selfridge_method_a_star(n: &BigInt) -> (BigInt, BigInt, BigInt) {
    let (d, p, q) = selfridge_method_a(n);
    if p.is_zero() {
        return (d, BigInt::zero(), BigInt::zero());
    }
    if q == BigInt::from(-1) {
        return (d, BigInt::from(5), BigInt::from(5));
    }
    (d, p, q)
}

/// Extra Strong parameter selection: Q=1, find P from 3 until Jacobi(P^2-4,n) = -1.
pub fn extra_strong_params(n: &BigInt) -> (BigInt, BigInt, BigInt) {
    let mut p = BigInt::from(3);
    loop {
        let d = p.clone() * p.clone() - BigInt::from(4);
        let j = jacobi_symbol(&d, n);
        if j == -1 {
            return (d, p, BigInt::one());
        } else if j == 0 {
            let g = d.abs().gcd(n);
            if g > BigInt::one() && g < *n {
                return (d, BigInt::zero(), BigInt::zero());
            }
            if d.abs() < *n {
                return (d, BigInt::zero(), BigInt::zero());
            }
        }
        p += 1;
        if p > BigInt::from(100000) {
            panic!("Extra Strong: no suitable P found within limit");
        }
    }
}

// ---------------------------------------------------------------------------
// Lucas probable prime tests
// ---------------------------------------------------------------------------

/// Regular Lucas PRP: U_{n+1} ≡ 0 (mod n).
pub fn is_regular_lucas_prp(n: &BigInt, _d: &BigInt, p: &BigInt, q: &BigInt) -> bool {
    let (u, _, _) = lucas_sequence(n, p, q, &(n + 1));
    u.is_zero()
}

/// Strong Lucas PRP.
pub fn is_strong_lucas_prp(n: &BigInt, _d: &BigInt, p: &BigInt, q: &BigInt) -> bool {
    let mut dd: BigInt = n + 1;
    let mut s = 0u32;
    while dd.is_even() {
        dd = dd / 2;
        s += 1;
    }
    let (u, v, qd) = lucas_sequence(n, p, q, &dd);
    if u.is_zero() {
        return true;
    }

    let mut vk = v;
    let mut qk = qd;
    for _ in 0..s {
        if vk.is_zero() {
            return true;
        }
        vk = mod_n(&(vk.clone() * vk.clone() - 2 * qk.clone()), n);
        qk = mod_n(&(qk.clone() * qk.clone()), n);
    }
    false
}

/// Extra Strong Lucas PRP (Q = 1).
pub fn is_extra_strong_lucas_prp(n: &BigInt, _d: &BigInt, p: &BigInt, q: &BigInt) -> bool {
    let mut dd: BigInt = n + 1;
    let mut s = 0u32;
    while dd.is_even() {
        dd = dd / 2;
        s += 1;
    }
    let (u, v, _) = lucas_sequence(n, p, q, &dd);
    if u.is_zero() && (v == BigInt::from(2) || v == n - 2) {
        return true;
    }

    let mut vk = v;
    for _ in 1..s {
        if vk.is_zero() {
            return true;
        }
        vk = mod_n(&(vk.clone() * vk.clone() - 2), n); // Q = 1
    }
    false
}

/// Almost Extra Strong Lucas PRP (Q = 1).
pub fn is_almost_extra_strong_lucas_prp(n: &BigInt, _d: &BigInt, p: &BigInt, q: &BigInt) -> bool {
    let mut dd: BigInt = n + 1;
    let mut s = 0u32;
    while dd.is_even() {
        dd = dd / 2;
        s += 1;
    }
    let (u, v, _) = lucas_sequence(n, p, q, &dd);
    if v == BigInt::from(2) || v == n - 2 {
        if u.is_zero() {
            return true;
        }
    }

    let mut vk = v;
    for _ in 0..s.saturating_sub(1) {
        if vk.is_zero() {
            return true;
        }
        vk = mod_n(&(vk.clone() * vk.clone() - 2), n); // Q = 1
    }
    false
}

// ---------------------------------------------------------------------------
// Strong probable prime with arbitrary base
// ---------------------------------------------------------------------------

fn strong_probable_prime_base(n: &BigInt, base: &BigInt) -> bool {
    if n < &BigInt::from(2) {
        return false;
    }
    if n.is_even() {
        return n == &BigInt::from(2);
    }

    let mut d: BigInt = n - 1;
    let mut s = 0u32;
    while d.is_even() {
        d = d / 2;
        s += 1;
    }

    let x = modpow(base, &d, n);
    if x == BigInt::one() || x == n - 1 {
        return true;
    }

    let mut x = x;
    for _ in 0..s - 1 {
        x = modpow(&x, &BigInt::from(2), n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Quick precheck
// ---------------------------------------------------------------------------

fn precheck(n: &BigInt) -> Option<bool> {
    if n < &BigInt::from(2) {
        return Some(false);
    }
    if n == &BigInt::from(2) {
        return Some(true);
    }
    if n.is_even() {
        return Some(false);
    }
    if n == &BigInt::from(3) {
        return Some(true);
    }
    if n % 3 == BigInt::zero() {
        return Some(false);
    }
    if is_square(n) {
        return Some(false);
    }
    None
}

// ---------------------------------------------------------------------------
// BPSW variants
// ---------------------------------------------------------------------------

/// Variant 2: Standard Strong BPSW — MR base 2 + Strong Lucas-Selfridge (Method A).
pub fn bpsw_standard(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let (d, p, q) = selfridge_method_a(n);
    if p.is_zero() {
        return false;
    }
    is_strong_lucas_prp(n, &d, &p, &q)
}

/// Variant 1: Regular Lucas BPSW — MR base 2 + Regular Lucas.
pub fn bpsw_regular_lucas(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let (d, p, q) = selfridge_method_a(n);
    if p.is_zero() {
        return false;
    }
    is_regular_lucas_prp(n, &d, &p, &q)
}

/// Variant 3: Extra Strong Lucas BPSW — MR base 2 + Extra Strong Lucas.
pub fn bpsw_extra_strong(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let (d, p, q) = extra_strong_params(n);
    if p.is_zero() {
        return false;
    }
    is_extra_strong_lucas_prp(n, &d, &p, &q)
}

/// Variant 4: Almost Extra Strong Lucas BPSW — MR base 2 + Almost Extra Strong.
pub fn bpsw_almost_extra_strong(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let (d, p, q) = extra_strong_params(n);
    if p.is_zero() {
        return false;
    }
    is_almost_extra_strong_lucas_prp(n, &d, &p, &q)
}

/// Variant 5: BFW / Enhanced BPSW (2021).
pub fn bfw_enhanced_bpsw(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let (d, p, q) = selfridge_method_a_star(n);
    if p.is_zero() {
        return false;
    }
    if !is_strong_lucas_prp(n, &d, &p, &q) {
        return false;
    }
    // Step 4: Lucas-V: V_{n+1} ≡ 2*Q (mod n)
    let (_, v, _) = lucas_sequence(n, &p, &q, &(n + 1));
    if v != mod_n(&(BigInt::from(2) * q.clone()), n) {
        return false;
    }
    // Step 5: Euler criterion: Q^((n+1)/2) ≡ Q * (Q/n) (mod n)
    let lhs = modpow(&q, &((n + 1) / 2), n);
    let j = jacobi_symbol(&q, n);
    let rhs = mod_n(&(q.clone() * BigInt::from(j)), n);
    if lhs != rhs {
        return false;
    }
    true
}

/// Variant 6: BPSW21 — same Lucas checks as BFW 2021.
pub fn bpsw21(n: &BigInt) -> bool {
    bfw_enhanced_bpsw(n)
}

/// Variant 7: Accelerated BFW — same conditions as BFW.
pub fn accelerated_bfw(n: &BigInt) -> bool {
    bfw_enhanced_bpsw(n)
}

/// Variant 8: SuperBFPSW (Hamburg 2025) — approximation (paper behind Cloudflare).
pub fn superbfpsw(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let (d, p, q) = selfridge_method_a(n);
    if p.is_zero() {
        return false;
    }
    // T derived from D (approximation; exact rule unknown without paper)
    let t: BigInt = d.abs() % n;
    let base = if t > BigInt::from(1) { t } else { BigInt::from(5) };
    if !strong_probable_prime_base(n, &base) {
        return false;
    }
    is_strong_lucas_prp(n, &d, &p, &q)
}

/// Variant 9: Fibonacci-BPSW — MR base 2 + Fibonacci when n ≡ 2,3 mod 5.
pub fn fibonacci_bpsw(n: &BigInt) -> bool {
    if let Some(r) = precheck(n) {
        return r;
    }
    if !is_strong_probable_prime_base2(n) {
        return false;
    }
    let n_mod_5 = n % 5;
    if n_mod_5 == BigInt::from(2) || n_mod_5 == BigInt::from(3) {
        // Fibonacci test: F_{n+1} ≡ 0 (mod n), P=1, Q=-1
        let (u, _, _) = lucas_sequence(n, &BigInt::from(1), &BigInt::from(-1), &(n + 1));
        u.is_zero()
    } else {
        let (d, p, q) = selfridge_method_a(n);
        if p.is_zero() {
            return false;
        }
        is_regular_lucas_prp(n, &d, &p, &q)
    }
}

// ---------------------------------------------------------------------------
// Variant registry
// ---------------------------------------------------------------------------

pub type VariantFunc = fn(&BigInt) -> bool;

pub const VARIANTS: [(u32, &'static str, VariantFunc); 9] = [
    (1, "Regular Lucas BPSW", bpsw_regular_lucas),
    (2, "Standard Strong BPSW", bpsw_standard),
    (3, "Extra Strong Lucas BPSW", bpsw_extra_strong),
    (4, "Almost Extra Strong Lucas BPSW", bpsw_almost_extra_strong),
    (5, "BFW / Enhanced BPSW (2021)", bfw_enhanced_bpsw),
    (6, "BPSW21", bpsw21),
    (7, "Accelerated BFW", accelerated_bfw),
    (8, "SuperBFPSW (Hamburg 2025)", superbfpsw),
    (9, "Fibonacci-BPSW", fibonacci_bpsw),
];

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn small_primes() -> Vec<BigInt> {
        vec![
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
            31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
            73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
        ]
        .into_iter()
        .map(BigInt::from)
        .collect()
    }

    fn small_composites() -> Vec<BigInt> {
        vec![
            4, 6, 8, 9, 10, 12, 14, 15, 16, 18,
            20, 21, 22, 24, 25, 26, 27, 28, 30, 32,
            33, 34, 35, 36, 38, 39, 40, 42, 44, 45,
            46, 48, 49, 50, 51, 52, 54, 55, 56, 57,
            58, 60, 62, 63, 64, 65, 66, 68, 69, 70,
            72, 74, 75, 76, 77, 78, 80, 81, 82, 84,
            85, 86, 87, 88, 90, 91, 92, 93, 94, 95,
            96, 98, 99, 100,
        ]
        .into_iter()
        .map(BigInt::from)
        .collect()
    }

    fn strong_psp_base2() -> Vec<BigInt> {
        vec![2047, 3277, 4033, 4681, 8321, 15841, 29341]
            .into_iter()
            .map(BigInt::from)
            .collect()
    }

    #[test]
    fn test_jacobi_symbol() {
        assert_eq!(jacobi_symbol(&BigInt::from(5), &BigInt::from(13)), -1);
        assert_eq!(jacobi_symbol(&BigInt::from(-7), &BigInt::from(13)), -1);
        assert_eq!(jacobi_symbol(&BigInt::from(1), &BigInt::from(7)), 1);
        assert_eq!(jacobi_symbol(&BigInt::from(0), &BigInt::from(7)), 0);
        assert_eq!(jacobi_symbol(&BigInt::from(2), &BigInt::from(7)), 1);
        assert_eq!(jacobi_symbol(&BigInt::from(2), &BigInt::from(5)), -1);
        assert_eq!(jacobi_symbol(&BigInt::from(3), &BigInt::from(11)), 1);
        assert_eq!(jacobi_symbol(&BigInt::from(4), &BigInt::from(9)), 1);
        assert_eq!(jacobi_symbol(&BigInt::from(6), &BigInt::from(9)), 0);
    }

    #[test]
    fn test_lucas_sequence_fibonacci() {
        let n = BigInt::from(1001);
        let (u, v, qk) = lucas_sequence(&n, &BigInt::from(1), &BigInt::from(-1), &BigInt::from(10));
        assert_eq!(u, BigInt::from(55));
        assert_eq!(v, BigInt::from(123));
        assert_eq!(qk, BigInt::from(1));
    }

    #[test]
    fn test_selfridge_method_a() {
        for p in small_primes().into_iter().skip(2) {
            let (d, pp, q) = selfridge_method_a(&p);
            assert_eq!(pp, BigInt::one());
            assert_eq!(jacobi_symbol(&d, &p), -1);
            assert_eq!(q, (BigInt::one() - d) / 4);
        }
    }

    #[test]
    fn test_selfridge_method_a_star() {
        for p in small_primes().into_iter().skip(2) {
            let (_, _, q) = selfridge_method_a_star(&p);
            assert_ne!(q, BigInt::from(-1));
        }
    }

    #[test]
    fn test_extra_strong_params() {
        for p in small_primes().into_iter().skip(2) {
            let (d, pp, q) = extra_strong_params(&p);
            assert_eq!(q, BigInt::one());
            assert!(pp > BigInt::from(2));
            assert_eq!(jacobi_symbol(&d, &p), -1);
            assert_eq!(d, pp.clone() * pp.clone() - BigInt::from(4));
        }
    }

    #[test]
    fn test_all_variants_accept_primes() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        for p in small_primes() {
            for f in &variants {
                assert!(f(&p), "variant rejected prime {}", p);
            }
        }
    }

    #[test]
    fn test_all_variants_reject_composites() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        for c in small_composites() {
            for f in &variants {
                assert!(!f(&c), "variant accepted composite {}", c);
            }
        }
    }

    #[test]
    fn test_strong_psp_base2_rejected() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        for n in strong_psp_base2() {
            for f in &variants {
                assert!(!f(&n), "variant accepted strong PSP base 2: {}", n);
            }
        }
    }

    #[test]
    fn test_edge_cases() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        for f in &variants {
            assert!(!f(&BigInt::zero()));
            assert!(!f(&BigInt::from(-5)));
            assert!(f(&BigInt::from(2)));
            assert!(!f(&BigInt::from(4)));
            assert!(!f(&BigInt::from(9)));
        }
    }

    #[test]
    fn test_variant_agreement_small_range() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        for n in 2..2001u32 {
            let bn = BigInt::from(n);
            let results: Vec<bool> = variants.iter().map(|f| f(&bn)).collect();
            let first = results[0];
            for r in &results[1..] {
                assert_eq!(*r, first, "variants disagree on n={}", n);
            }
        }
    }

    #[test]
    fn test_large_primes() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        let primes = vec![
            BigInt::from(1000003),
            BigInt::from(1000033),
            BigInt::from(1000000007),
            BigInt::from(1000000009),
        ];
        for p in primes {
            for f in &variants {
                assert!(f(&p), "variant rejected prime {}", p);
            }
        }
    }

    #[test]
    fn test_large_composites() {
        let variants: Vec<VariantFunc> = VARIANTS.iter().map(|(_, _, f)| *f).collect();
        let composites = vec![
            BigInt::from(1000001),
            BigInt::from(1000000001),
        ];
        for c in composites {
            for f in &variants {
                assert!(!f(&c), "variant accepted composite {}", c);
            }
        }
    }
}
