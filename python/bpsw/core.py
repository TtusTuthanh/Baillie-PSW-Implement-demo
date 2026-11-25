"""Core mathematical primitives for BPSW variants.

Implements:
  - Jacobi symbol (a/n)
  - Perfect square detection
  - Strong probable prime test base 2 (Miller-Rabin)
  - Lucas sequence U_k, V_k, Q^k via binary ladder
  - Selfridge parameter selection: Method A, Method A*, Extra Strong (Q=1)
"""

from __future__ import annotations

import math


def jacobi_symbol(a: int, n: int) -> int:
    """Compute the Jacobi symbol (a/n). n must be odd and positive.

    Handles negative a by reducing modulo n first.
    Returns 0 if gcd(a, n) > 1, +1 or -1 otherwise.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be odd and positive")
    a = a % n
    if a == 0:
        return 0
    result = 1
    m = n  # working copy of modulus
    while a:
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5):
                result = -result
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3:
            result = -result
        a %= m
    return result if m == 1 else 0


def is_square(n: int) -> bool:
    """Check if n is a perfect square."""
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def is_strong_probable_prime_base2(n: int) -> bool:
    """Strong probable prime test base 2 (Miller-Rabin).

    Returns True if n is a strong probable prime to base 2,
    False if n is definitely composite.
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Write n - 1 = d * 2^s with d odd
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(2, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return True
    return False


def lucas_sequence(n: int, P: int, Q: int, k: int) -> tuple[int, int, int]:
    """Compute (U_k mod n, V_k mod n, Q^k mod n) using the binary ladder.

    Lucas sequences:
      U_0 = 0,  U_1 = 1,  U_k = P * U_{k-1} - Q * U_{k-2}
      V_0 = 2,  V_1 = P,  V_k = P * V_{k-1} - Q * V_{k-2}

    n: odd positive modulus
    P, Q: Lucas parameters (integers, may be negative)
    k: non-negative index

    Returns (U_k, V_k, Q^k) all reduced mod n.
    """
    if k == 0:
        return (0 % n, 2 % n, 1 % n)

    inv2 = (n + 1) // 2  # modular inverse of 2 mod n (n is odd)
    D = (P * P - 4 * Q) % n  # discriminant

    # Binary representation of k, MSB first
    bits = []
    kk = k
    while kk:
        bits.append(kk & 1)
        kk >>= 1
    bits.reverse()

    # State for current index m: (U_m, V_m, Q^m)
    U = 0   # U_0
    V = 2   # V_0
    Qk = 1  # Q^0

    for bit in bits:
        # Double: m -> 2m
        #   U_{2m} = U_m * V_m
        #   V_{2m} = V_m^2 - 2 * Q^m
        #   Q^{2m} = (Q^m)^2
        U = (U * V) % n
        V = (V * V - 2 * Qk) % n
        Qk = (Qk * Qk) % n

        if bit:
            # Advance: 2m -> 2m + 1
            #   U_{2m+1} = (P * U_{2m} + V_{2m}) / 2   mod n
            #   V_{2m+1} = (D * U_{2m} + P * V_{2m}) / 2   mod n
            #   Q^{2m+1} = Q^{2m} * Q   mod n
            U_new = (P * U + V) * inv2 % n
            V_new = (D * U + P * V) * inv2 % n
            Qk_new = (Qk * Q) % n
            U, V, Qk = U_new, V_new, Qk_new

    return (U, V, Qk)


def selfridge_method_a(n: int) -> tuple[int, int, int]:
    """Selfridge Method A parameter selection.

    Scan D through 5, -7, 9, -11, 13, -15, ... until Jacobi(D, n) == -1.
    Set P = 1, Q = (1 - D) / 4.

    Returns (D, P, Q).
    Returns (D, 0, 0) as a composite signal when a nontrivial factor is found.
    """
    D = 5
    for _ in range(100000):  # safety limit
        j = jacobi_symbol(D, n)
        if j == -1:
            P = 1
            Q = (1 - D) // 4
            return (D, P, Q)
        elif j == 0:
            # D shares a factor with n
            g = math.gcd(abs(D), n)
            if 1 < g < n:
                return (D, 0, 0)  # composite: found nontrivial factor
            if abs(D) < n:
                return (D, 0, 0)  # composite: D < n and shares factor
        # Next D in sequence: 5, -7, 9, -11, 13, -15, ...
        if D > 0:
            D = -(D + 2)
        else:
            D = -D + 2
    raise ValueError("Selfridge Method A: no suitable D found within limit")


def selfridge_method_a_star(n: int) -> tuple[int, int, int]:
    """Selfridge Method A* (BFW 2021 refinement).

    Same as Method A, but if Q == -1 (which happens when D == 5),
    set P = Q = 5 instead. This forces Q not congruent to ±1 (mod n),
    which dramatically reduces simultaneous Lucas/V pseudoprimes.

    Returns (D, P, Q).
    Returns (D, 0, 0) as a composite signal.
    """
    D, P, Q = selfridge_method_a(n)
    if P == 0:
        return (D, 0, 0)
    if Q == -1:
        P = 5
        Q = 5
    return (D, P, Q)


def extra_strong_params(n: int) -> tuple[int, int, int]:
    """Extra Strong Lucas parameter selection (Q = 1).

    Find P starting from 3 such that Jacobi(P^2 - 4, n) == -1.
    Set D = P^2 - 4, Q = 1.

    Returns (D, P, Q).
    Returns (D, 0, 0) as a composite signal.
    """
    P = 3
    while True:
        D = P * P - 4
        j = jacobi_symbol(D, n)
        if j == -1:
            return (D, P, 1)
        elif j == 0:
            g = math.gcd(abs(D), n)
            if 1 < g < n:
                return (D, 0, 0)
            if abs(D) < n:
                return (D, 0, 0)
        P += 1
        if P > 100000:
            raise ValueError("Extra Strong: no suitable P found within limit")
