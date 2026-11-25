"""Tests for BPSW variants.

Covers:
  - Small primes (all variants must accept)
  - Small composites (all variants must reject)
  - Strong pseudoprimes base 2 (must be caught by Lucas test)
  - Jacobi symbol correctness
  - Lucas sequence correctness (Fibonacci numbers)
  - Selfridge parameter selection correctness
  - Edge cases (0, 1, negative, even, squares)
  - Variant agreement on small range
  - Known primes up to 10^18
"""

from __future__ import annotations

import pytest

from bpsw import (
    bpsw_standard,
    bpsw_regular_lucas,
    bpsw_extra_strong,
    bpsw_almost_extra_strong,
    is_bfw_enhanced_bpsw,
    is_bpsw21,
    is_accelerated_bfw,
    is_superbfpsw,
    is_fibonacci_bpsw,
    VARIANTS,
    jacobi_symbol,
    is_square,
    is_strong_probable_prime_base2,
    lucas_sequence,
    is_regular_lucas_prp,
    is_strong_lucas_prp,
    is_extra_strong_lucas_prp,
    is_almost_extra_strong_lucas_prp,
    selfridge_method_a,
    selfridge_method_a_star,
    extra_strong_params,
)

# All variant functions (for iteration)
ALL_VARIANTS = [func for _, func in VARIANTS.values()]

# ---------------------------------------------------------------------------
# Known small primes
# ---------------------------------------------------------------------------

SMALL_PRIMES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199,
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
    269, 271, 277, 281, 283, 293,
]

# ---------------------------------------------------------------------------
# Known small composites
# ---------------------------------------------------------------------------

SMALL_COMPOSITES = [
    4, 6, 8, 9, 10, 12, 14, 15, 16, 18,
    20, 21, 22, 24, 25, 26, 27, 28, 30, 32,
    33, 34, 35, 36, 38, 39, 40, 42, 44, 45,
    46, 48, 49, 50, 51, 52, 54, 55, 56, 57,
    58, 60, 62, 63, 64, 65, 66, 68, 69, 70,
    72, 74, 75, 76, 77, 78, 80, 81, 82, 84,
    85, 86, 87, 88, 90, 91, 92, 93, 94, 95,
    96, 98, 99, 100,
    121, 143, 169, 187, 209, 221, 247, 253, 289, 299,
    323, 361, 377, 391, 437, 493, 527, 529, 551, 589,
]

# ---------------------------------------------------------------------------
# Strong pseudoprimes base 2 (pass MR base 2, should fail Lucas)
# ---------------------------------------------------------------------------

STRONG_PSP_BASE2 = [2047, 3277, 4033, 4681, 8321, 15841, 29341]

# ---------------------------------------------------------------------------
# Larger primes (up to ~10^18)
# ---------------------------------------------------------------------------

LARGE_PRIMES = [
    1000003,
    1000033,
    1000037,
    1000000007,
    1000000009,
    1000000000000000003,
    1000000000000000009,
]

# ---------------------------------------------------------------------------
# Larger composites
# ---------------------------------------------------------------------------

LARGE_COMPOSITES = [
    1000001,   # 101 * 9901
    1000000001,  # 7 * 11 * 13 * 19 * 52579
    1000000000000000001,  # composite
]


# ===========================================================================
# Basic correctness: all variants accept primes, reject composites
# ===========================================================================

@pytest.mark.parametrize("p", SMALL_PRIMES)
def test_all_variants_accept_primes(p):
    """All BPSW variants should accept small primes."""
    for name, func in VARIANTS.values():
        assert func(p) is True, f"{name} rejected prime {p}"


@pytest.mark.parametrize("c", SMALL_COMPOSITES)
def test_all_variants_reject_composites(c):
    """All BPSW variants should reject small composites."""
    for name, func in VARIANTS.values():
        assert func(c) is False, f"{name} accepted composite {c}"


@pytest.mark.parametrize("p", LARGE_PRIMES)
def test_all_variants_accept_large_primes(p):
    """All BPSW variants should accept larger primes."""
    for name, func in VARIANTS.values():
        assert func(p) is True, f"{name} rejected prime {p}"


@pytest.mark.parametrize("c", LARGE_COMPOSITES)
def test_all_variants_reject_large_composites(c):
    """All BPSW variants should reject larger composites."""
    for name, func in VARIANTS.values():
        assert func(c) is False, f"{name} accepted composite {c}"


# ===========================================================================
# Strong PSP base 2 should be caught by all variants
# ===========================================================================

@pytest.mark.parametrize("n", STRONG_PSP_BASE2)
def test_strong_psp_base2_rejected(n):
    """Strong pseudoprimes base 2 should be rejected by all BPSW variants."""
    for name, func in VARIANTS.values():
        assert func(n) is False, f"{name} accepted strong PSP base 2: {n}"


# ===========================================================================
# Jacobi symbol tests
# ===========================================================================

def test_jacobi_symbol():
    """Test Jacobi symbol against known values."""
    assert jacobi_symbol(5, 13) == -1
    assert jacobi_symbol(-7, 13) == -1
    assert jacobi_symbol(1, 7) == 1
    assert jacobi_symbol(0, 7) == 0
    assert jacobi_symbol(2, 7) == 1   # 7 ≡ -1 mod 8, so (2/7) = 1
    assert jacobi_symbol(2, 5) == -1  # 5 ≡ -3 mod 8, so (2/5) = -1
    assert jacobi_symbol(3, 11) == 1
    # (4/9) = 1 because 4 = 2^2 is a quadratic residue mod 9
    assert jacobi_symbol(4, 9) == 1
    # (6/9) = 0 because gcd(6,9) = 3 > 1
    assert jacobi_symbol(6, 9) == 0


# ===========================================================================
# Lucas sequence tests
# ===========================================================================

def test_lucas_sequence_fibonacci():
    """Test Lucas sequence with P=1, Q=-1 gives Fibonacci/Lucas numbers.

    F: 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, ...
    L: 2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, ...

    Uses odd modulus n=1001 so that inv2 exists.
    """
    U, V, Qk = lucas_sequence(1001, 1, -1, 10)
    assert U == 55   # F_10 = 55
    assert V == 123  # L_10 = 123
    assert Qk == 1   # (-1)^10 = 1

    U, V, Qk = lucas_sequence(1001, 1, -1, 0)
    assert U == 0
    assert V == 2
    assert Qk == 1


def test_lucas_sequence_identity():
    """Verify Lucas sequence identities for various parameters.

    Uses odd modulus n=10001 so inv2 exists.
    P=3, Q=2: D = 9-8 = 1
    U: 0, 1, 3, 7, 15, 31, ... (U_k = 2^k - 1)
    V: 2, 3, 5, 9, 17, 33, ... (V_k = 2^k + 1)
    """
    for k in range(0, 8):
        U, V, Qk = lucas_sequence(10001, 3, 2, k)
        assert U == (2 ** k - 1) % 10001, f"U_{k} mismatch"
        assert V == (2 ** k + 1) % 10001, f"V_{k} mismatch"
        assert Qk == pow(2, k, 10001), f"Q^{k} mismatch"


# ===========================================================================
# Selfridge parameter selection
# ===========================================================================

def test_selfridge_method_a():
    """Test Selfridge Method A: Jacobi(D, n) == -1 for selected D."""
    for p in SMALL_PRIMES[2:]:  # skip 2, 3
        D, P, Q = selfridge_method_a(p)
        assert P == 1
        assert jacobi_symbol(D, p) == -1
        assert Q == (1 - D) // 4


def test_selfridge_method_a_star():
    """Test Selfridge Method A*: Q != -1 always."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = selfridge_method_a_star(p)
        assert jacobi_symbol(D, p) == -1
        assert Q != -1, f"Method A* produced Q=-1 for n={p}"


def test_extra_strong_params():
    """Test Extra Strong parameter selection: Q=1, Jacobi(D,n) == -1."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = extra_strong_params(p)
        assert Q == 1
        assert P >= 3
        assert jacobi_symbol(D, p) == -1
        assert D == P * P - 4


# ===========================================================================
# Edge cases
# ===========================================================================

@pytest.mark.parametrize("func", ALL_VARIANTS)
def test_edge_cases(func):
    """All variants must handle edge cases correctly."""
    assert func(0) is False
    assert func(1) is False
    assert func(-5) is False
    assert func(2) is True
    assert func(4) is False   # even composite
    assert func(9) is False   # square
    assert func(25) is False  # square
    assert func(49) is False  # square


# ===========================================================================
# Variant agreement on range 2..2000
# ===========================================================================

def test_variant_agreement_small_range():
    """All variants should agree on primality for 2..2000.
    
    Note: Fibonacci-BPSW uses a different Lucas test for n≡2,3 mod 5,
    but for small numbers all variants should agree since the Lucas test
    is strong enough to catch all composites in this range.
    """
    results_by_n = {}
    for n in range(2, 2001):
        results = [func(n) for _, func in VARIANTS.values()]
        results_by_n[n] = results
        # All variants must agree
        first = results[0]
        for r in results[1:]:
            assert r == first, f"Variants disagree on n={n}: {results}"

    # Cross-check with simple trial division for the small range
    def is_prime_trial(n):
        if n < 2:
            return False
        if n % 2 == 0:
            return n == 2
        i = 3
        while i * i <= n:
            if n % i == 0:
                return False
            i += 2
        return True

    for n in range(2, 2001):
        expected = is_prime_trial(n)
        got = results_by_n[n][0]
        assert got == expected, f"Variant result {got} != trial division {expected} for n={n}"


# ===========================================================================
# BFW-specific checks
# ===========================================================================

def test_bfw_method_a_star_q_not_minus_one():
    """BFW uses Method A* which ensures Q != -1."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = selfridge_method_a_star(p)
        assert Q != -1, f"Method A* Q=-1 for n={p}"


def test_bfw_euler_criterion():
    """Verify Euler criterion condition on primes."""
    for p in SMALL_PRIMES[5:]:  # skip small ones
        D, P, Q = selfridge_method_a_star(p)
        if P == 0:
            continue
        lhs = pow(Q, (p + 1) // 2, p)
        rhs = (Q * jacobi_symbol(Q, p)) % p
        assert lhs == rhs, f"Euler criterion failed for p={p}, Q={Q}"


# ===========================================================================
# Fibonacci-BPSW specific
# ===========================================================================

@pytest.mark.parametrize("p", SMALL_PRIMES)
def test_fibonacci_bpsw_primes(p):
    """Fibonacci-BPSW accepts all small primes."""
    assert is_fibonacci_bpsw(p) is True


@pytest.mark.parametrize("c", SMALL_COMPOSITES)
def test_fibonacci_bpsw_composites(c):
    """Fibonacci-BPSW rejects all small composites."""
    assert is_fibonacci_bpsw(c) is False


def test_fibonacci_bpsw_n_mod_5():
    """For n ≡ 2,3 mod 5, Fibonacci-BPSW uses Fibonacci test."""
    # Just verify it runs without error for a range
    for n in range(2, 500):
        is_fibonacci_bpsw(n)  # should not raise


# ===========================================================================
# Lucas PRP individual tests
# ===========================================================================

def test_strong_lucas_prp_known_prime():
    """Strong Lucas PRP accepts known primes with Selfridge params."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = selfridge_method_a(p)
        assert is_strong_lucas_prp(p, D, P, Q), f"Strong Lucas rejected prime {p}"


def test_extra_strong_lucas_prp_known_prime():
    """Extra Strong Lucas PRP accepts known primes."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = extra_strong_params(p)
        assert is_extra_strong_lucas_prp(p, D, P, Q), f"Extra Strong Lucas rejected prime {p}"


def test_almost_extra_strong_lucas_prp_known_prime():
    """Almost Extra Strong Lucas PRP accepts known primes."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = extra_strong_params(p)
        assert is_almost_extra_strong_lucas_prp(p, D, P, Q), f"Almost Extra Strong Lucas rejected prime {p}"


def test_regular_lucas_prp_known_prime():
    """Regular Lucas PRP accepts known primes."""
    for p in SMALL_PRIMES[2:]:
        D, P, Q = selfridge_method_a(p)
        assert is_regular_lucas_prp(p, D, P, Q), f"Regular Lucas rejected prime {p}"
