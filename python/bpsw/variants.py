"""BPSW variants: all documented configurations.

Each variant returns True (probable prime) or False (composite).
All variants reject even numbers, n < 2, and perfect squares up front.
"""

from __future__ import annotations

from .core import (
    jacobi_symbol,
    is_square,
    is_strong_probable_prime_base2,
    lucas_sequence,
    selfridge_method_a,
    selfridge_method_a_star,
    extra_strong_params,
)


# ---------------------------------------------------------------------------
# Quick precheck
# ---------------------------------------------------------------------------

def _precheck(n: int) -> bool | None:
    """Quick checks: n<2, even, small primes, perfect squares.

    Returns True if definitely prime, False if definitely composite,
    None if the full test is needed.
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    if n == 3:
        return True
    if n % 3 == 0:
        return False
    if is_square(n):
        return False
    return None


# ---------------------------------------------------------------------------
# Individual Lucas probable prime tests
# ---------------------------------------------------------------------------

def is_regular_lucas_prp(n: int, D: int, P: int, Q: int) -> bool:
    """Regular Lucas probable prime test: U_{n+1} ≡ 0 (mod n).

    Weaker than Strong Lucas; used in "Regular Lucas BPSW" variant.
    """
    U, _, _ = lucas_sequence(n, P, Q, n + 1)
    return U == 0


def is_strong_lucas_prp(n: int, D: int, P: int, Q: int) -> bool:
    """Strong Lucas probable prime test.

    Write n + 1 = d * 2^s with d odd. Accept if:
      U_d ≡ 0 (mod n),  OR
      V_{d·2^r} ≡ 0 (mod n) for some 0 ≤ r < s.
    """
    d = n + 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    U, V, Qd = lucas_sequence(n, P, Q, d)
    if U == 0:
        return True
    Vk = V
    Qk = Qd
    for _ in range(s):
        if Vk == 0:
            return True
        Vk = (Vk * Vk - 2 * Qk) % n
        Qk = (Qk * Qk) % n
    return False


def is_extra_strong_lucas_prp(n: int, D: int, P: int, Q: int) -> bool:
    """Extra Strong Lucas probable prime test (Q = 1).

    Write n + 1 = d * 2^s with d odd. Accept if:
      U_d ≡ 0 AND V_d ≡ ±2 (mod n),  OR
      V_{d·2^r} ≡ 0 (mod n) for some 1 ≤ r < s.
    """
    d = n + 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    U, V, _ = lucas_sequence(n, P, Q, d)
    if U == 0 and (V == 2 or V == n - 2):
        return True
    Vk = V
    # Q = 1, so Q^k = 1 for all k
    for _ in range(1, s):
        if Vk == 0:
            return True
        Vk = (Vk * Vk - 2) % n
    return False


def is_almost_extra_strong_lucas_prp(n: int, D: int, P: int, Q: int) -> bool:
    """Almost Extra Strong Lucas probable prime test (Q = 1).

    Like Extra Strong but skips the standalone U_d ≡ 0 check.
    Accept if:
      V_d ≡ ±2 AND U_d ≡ 0 (mod n),  OR
      V_{d·2^r} ≡ 0 (mod n) for some 0 ≤ r < s - 1.
    """
    d = n + 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    U, V, _ = lucas_sequence(n, P, Q, d)
    if V == 2 or V == n - 2:
        if U == 0:
            return True
    Vk = V
    # Q = 1, so Q^k = 1 for all k
    for _ in range(s - 1):
        if Vk == 0:
            return True
        Vk = (Vk * Vk - 2) % n
    return False


# ---------------------------------------------------------------------------
# Strong probable prime test with arbitrary base
# ---------------------------------------------------------------------------

def _strong_probable_prime_base(n: int, base: int) -> bool:
    """Strong probable prime test with given base."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(base % n, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return True
    return False


# ---------------------------------------------------------------------------
# BPSW variants
# ---------------------------------------------------------------------------

def bpsw_standard(n: int) -> bool:
    """Variant 2: Standard Strong BPSW.

    Strong Miller–Rabin base 2 + Strong Lucas-Selfridge (Method A).
    This is the most common modern BPSW configuration.
    """
    result = _precheck(n)
    if result is not None:
        return result
    if not is_strong_probable_prime_base2(n):
        return False
    D, P, Q = selfridge_method_a(n)
    if P == 0:
        return False
    return is_strong_lucas_prp(n, D, P, Q)


def bpsw_regular_lucas(n: int) -> bool:
    """Variant 1: Regular Lucas BPSW.

    Strong Miller–Rabin base 2 + Regular Lucas (U_{n+1} ≡ 0).
    Weaker than Strong Lucas BPSW.
    """
    result = _precheck(n)
    if result is not None:
        return result
    if not is_strong_probable_prime_base2(n):
        return False
    D, P, Q = selfridge_method_a(n)
    if P == 0:
        return False
    return is_regular_lucas_prp(n, D, P, Q)


def bpsw_extra_strong(n: int) -> bool:
    """Variant 3: Extra Strong Lucas BPSW.

    Strong Miller–Rabin base 2 + Extra Strong Lucas (Q=1, find P).
    """
    result = _precheck(n)
    if result is not None:
        return result
    if not is_strong_probable_prime_base2(n):
        return False
    D, P, Q = extra_strong_params(n)
    if P == 0:
        return False
    return is_extra_strong_lucas_prp(n, D, P, Q)


def bpsw_almost_extra_strong(n: int) -> bool:
    """Variant 4: Almost Extra Strong Lucas BPSW.

    Strong Miller–Rabin base 2 + Almost Extra Strong Lucas.
    Used by Go math/big (Baillie-OEIS method C).
    """
    result = _precheck(n)
    if result is not None:
        return result
    if not is_strong_probable_prime_base2(n):
        return False
    D, P, Q = extra_strong_params(n)
    if P == 0:
        return False
    return is_almost_extra_strong_lucas_prp(n, D, P, Q)


def is_bfw_enhanced_bpsw(n: int) -> bool:
    """Variant 5: BFW / Enhanced BPSW (Baillie-Fiori-Wagstaff 2021).

    5 steps:
      1. Strong base-2 probable prime test
      2. Method A* parameter selection
      3. Strong Lucas probable prime test
      4. Lucas-V condition: V_{n+1} ≡ 2·Q (mod n)
      5. Euler criterion: Q^((n+1)/2) ≡ Q · (Q/n) (mod n)

    All steps must pass.
    """
    result = _precheck(n)
    if result is not None:
        return result

    # Step 1: Strong base-2 MR
    if not is_strong_probable_prime_base2(n):
        return False

    # Step 2: Method A* parameter selection
    D, P, Q = selfridge_method_a_star(n)
    if P == 0:
        return False

    # Step 3: Strong Lucas
    if not is_strong_lucas_prp(n, D, P, Q):
        return False

    # Step 4: Lucas-V condition: V_{n+1} ≡ 2·Q (mod n)
    _, V, _ = lucas_sequence(n, P, Q, n + 1)
    if V != (2 * Q) % n:
        return False

    # Step 5: Euler criterion: Q^((n+1)/2) ≡ Q · (Q/n) (mod n)
    lhs = pow(Q, (n + 1) // 2, n)
    rhs = (Q * jacobi_symbol(Q, n)) % n
    if lhs != rhs:
        return False

    return True


def is_bpsw21(n: int) -> bool:
    """Variant 6: BPSW21.

    Same Lucas checks as BFW 2021 (Strong Lucas + Lucas-V + Euler criterion).
    In the crypto_primes library, Bpsw21 is a LucasCheck mode, not a
    standalone algorithm; here it is exposed as a full BPSW variant.
    """
    return is_bfw_enhanced_bpsw(n)


def is_accelerated_bfw(n: int) -> bool:
    """Variant 7: Accelerated BFW.

    Same mathematical acceptance conditions as BFW 2021. The "acceleration"
    refers to optimized Lucas sequence computation (same congruences,
    faster constant factors). Here it delegates to the BFW implementation.
    """
    return is_bfw_enhanced_bpsw(n)


def is_superbfpsw(n: int) -> bool:
    """Variant 8: SuperBFPSW (Mike Hamburg, 2025–2026).

    NOTE: The original paper (ePrint 2025/2083) was not accessible
    (Cloudflare block on ePrint). Only the abstract was obtained:
      "SuperBFPSW is stronger than Baillie-Fiori-Wagstaff, but also
       faster than the original Baillie-PSW."

    This implementation is an approximation based on the published
    description:
      - Strong Fermat base 2
      - Lucas parameter selection (Method A, approximation)
      - Strong Fermat base T (T derived from D; base 5 in special case)
      - Strong Lucas test

    Exact algorithm details (step 4 T determination, accelerated Lucas
    specifics, Euler-Lucas removal) require the paper. Results from this
    variant should be treated as approximate, not authoritative.
    """
    result = _precheck(n)
    if result is not None:
        return result

    # Step 1: Strong Fermat base 2
    if not is_strong_probable_prime_base2(n):
        return False

    # Step 2-3: Find discriminant D and parameter T
    # Using Method A as approximation (paper details unavailable)
    D, P, Q = selfridge_method_a(n)
    if P == 0:
        return False

    # T derived from D (approximation; exact rule unknown without paper)
    T = abs(D) % n
    if T <= 1:
        T = 5  # special case fallback: base 5

    # Step 4: Strong Fermat base T (or base 5 in special case)
    if not _strong_probable_prime_base(n, T):
        return False

    # Step 5: Accelerated Strong Lucas (same conditions, approximate)
    if not is_strong_lucas_prp(n, D, P, Q):
        return False

    return True


def is_fibonacci_bpsw(n: int) -> bool:
    """Variant 9: Fibonacci-BPSW (Chen & Greene 2003 / FLINT).

    Strong Miller–Rabin base 2 + Lucas/Fibonacci test:
      - If n ≡ 2 or 3 (mod 5): use Fibonacci test (P=1, Q=−1),
        check F_{n+1} ≡ 0 (mod n), i.e. U_{n+1} ≡ 0.
      - Otherwise: use regular Lucas-Selfridge test.
    """
    result = _precheck(n)
    if result is not None:
        return result

    # Step 1: Strong base-2 MR
    if not is_strong_probable_prime_base2(n):
        return False

    # Step 2: Lucas or Fibonacci depending on n mod 5
    if n % 5 in (2, 3):
        # Fibonacci test: F_{n+1} ≡ 0 mod n (P=1, Q=−1)
        U, _, _ = lucas_sequence(n, 1, -1, n + 1)
        return U == 0
    else:
        # Regular Lucas-Selfridge
        D, P, Q = selfridge_method_a(n)
        if P == 0:
            return False
        return is_regular_lucas_prp(n, D, P, Q)


# ---------------------------------------------------------------------------
# Variant registry
# ---------------------------------------------------------------------------

VARIANTS: dict[int, tuple[str, callable]] = {
    1: ("Regular Lucas BPSW", bpsw_regular_lucas),
    2: ("Standard Strong BPSW", bpsw_standard),
    3: ("Extra Strong Lucas BPSW", bpsw_extra_strong),
    4: ("Almost Extra Strong Lucas BPSW", bpsw_almost_extra_strong),
    5: ("BFW / Enhanced BPSW (2021)", is_bfw_enhanced_bpsw),
    6: ("BPSW21", is_bpsw21),
    7: ("Accelerated BFW", is_accelerated_bfw),
    8: ("SuperBFPSW (Hamburg 2025)", is_superbfpsw),
    9: ("Fibonacci-BPSW", is_fibonacci_bpsw),
}
