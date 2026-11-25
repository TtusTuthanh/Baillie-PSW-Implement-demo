#include "bpsw.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Core mathematical primitives
// ---------------------------------------------------------------------------

Int mod_n(Int x, Int n) {
    // Normalize x mod n to [0, n-1]. __int128 % can return negative.
    Int r = x % n;
    if (r < 0) r += n;
    return r;
}

int jacobi_symbol(Int a, Int n) {
    if (n <= 0 || n % 2 == 0)
        throw std::invalid_argument("n must be odd and positive");
    a = mod_n(a, n);
    if (a == 0) return 0;
    int result = 1;
    Int m = n;
    while (a != 0) {
        while (a % 2 == 0) {
            a /= 2;
            Int r = m % 8;
            if (r == 3 || r == 5) result = -result;
        }
        std::swap(a, m);
        if (a % 4 == 3 && m % 4 == 3) result = -result;
        a = a % m;
    }
    return (m == 1) ? result : 0;
}

bool is_square(Int n) {
    if (n < 0) return false;
    // Integer square root via Newton's method
    Int x = n;
    while (true) {
        Int y = (x + n / x) / 2;
        if (y >= x) break;
        x = y;
    }
    return x * x == n;
}

Int modpow(Int base, Int exp, Int modulus) {
    if (modulus == 1) return 0;
    base = mod_n(base, modulus);
    Int result = 1;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % modulus;
        base = (base * base) % modulus;
        exp >>= 1;
    }
    return result;
}

bool is_strong_probable_prime_base2(Int n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;

    // Write n - 1 = d * 2^s with d odd
    Int d = n - 1;
    int s = 0;
    while (d % 2 == 0) { d /= 2; s++; }

    Int x = modpow(2, d, n);
    if (x == 1 || x == n - 1) return true;

    Int cur = x;
    for (int i = 0; i < s - 1; i++) {
        cur = (cur * cur) % n;
        if (cur == n - 1) return true;
    }
    return false;
}

std::tuple<Int, Int, Int> lucas_sequence(Int n, Int P, Int Q, Int k) {
    if (k == 0) {
        return { mod_n(0, n), mod_n(2, n), mod_n(1, n) };
    }

    Int inv2 = (n + 1) / 2;   // modular inverse of 2 mod n (n is odd)
    Int D = mod_n(P * P - 4 * Q, n);

    // Binary representation of k, MSB first
    std::vector<int> bits;
    Int kk = k;
    while (kk > 0) {
        bits.push_back((int)(kk & 1));
        kk >>= 1;
    }
    // bits are LSB first; reverse to get MSB first
    std::reverse(bits.begin(), bits.end());

    Int U = 0;   // U_0
    Int V = 2;   // V_0
    Int Qk = 1;  // Q^0

    for (int bit : bits) {
        // Double: m -> 2m
        U = mod_n(U * V, n);
        V = mod_n(V * V - 2 * Qk, n);
        Qk = mod_n(Qk * Qk, n);

        if (bit == 1) {
            // Advance: 2m -> 2m+1
            Int U_new = mod_n((P * U + V) * inv2, n);
            Int V_new = mod_n((D * U + P * V) * inv2, n);
            Int Qk_new = mod_n(Qk * Q, n);
            U = U_new;
            V = V_new;
            Qk = Qk_new;
        }
    }

    return { U, V, Qk };
}

std::tuple<Int, Int, Int> selfridge_method_a(Int n) {
    Int D = 5;
    for (int i = 0; i < 100000; i++) {
        int j = jacobi_symbol(D, n);
        if (j == -1) {
            Int P = 1;
            Int Q = (1 - D) / 4;
            return { D, P, Q };
        } else if (j == 0) {
            // D shares a factor with n
            // Use Euclidean algorithm for gcd
            Int a = D < 0 ? -D : D;
            Int b = n;
            while (b != 0) { Int t = b; b = a % b; a = t; }
            Int g = a;
            if (g > 1 && g < n) return { D, 0, 0 };
            if (a < n) return { D, 0, 0 };
        }
        // Next D: 5, -7, 9, -11, ...
        if (D > 0) D = -(D + 2);
        else D = -D + 2;
    }
    throw std::runtime_error("Selfridge Method A: no suitable D found");
}

std::tuple<Int, Int, Int> selfridge_method_a_star(Int n) {
    auto [D, P, Q] = selfridge_method_a(n);
    if (P == 0) return { D, 0, 0 };
    if (Q == -1) { P = 5; Q = 5; }
    return { D, P, Q };
}

std::tuple<Int, Int, Int> extra_strong_params(Int n) {
    Int P = 3;
    while (true) {
        Int D = P * P - 4;
        int j = jacobi_symbol(D, n);
        if (j == -1) return { D, P, 1 };
        else if (j == 0) {
            Int a = D < 0 ? -D : D;
            Int b = n;
            while (b != 0) { Int t = b; b = a % b; a = t; }
            Int g = a;
            if (g > 1 && g < n) return { D, 0, 0 };
            if (a < n) return { D, 0, 0 };
        }
        P++;
        if (P > 100000) throw std::runtime_error("Extra Strong: no suitable P found");
    }
}

// ---------------------------------------------------------------------------
// Lucas probable prime tests
// ---------------------------------------------------------------------------

bool is_regular_lucas_prp(Int n, Int D, Int P, Int Q) {
    (void)D; // unused
    auto [U, _, __] = lucas_sequence(n, P, Q, n + 1);
    return U == 0;
}

bool is_strong_lucas_prp(Int n, Int D, Int P, Int Q) {
    (void)D; // unused
    Int dd = n + 1;
    int s = 0;
    while (dd % 2 == 0) { dd /= 2; s++; }
    auto [U, V, Qd] = lucas_sequence(n, P, Q, dd);
    if (U == 0) return true;

    Int Vk = V;
    Int Qk = Qd;
    for (int i = 0; i < s; i++) {
        if (Vk == 0) return true;
        Vk = mod_n(Vk * Vk - 2 * Qk, n);
        Qk = mod_n(Qk * Qk, n);
    }
    return false;
}

bool is_extra_strong_lucas_prp(Int n, Int D, Int P, Int Q) {
    (void)D; // unused
    Int dd = n + 1;
    int s = 0;
    while (dd % 2 == 0) { dd /= 2; s++; }
    auto [U, V, __] = lucas_sequence(n, P, Q, dd);
    if (U == 0 && (V == 2 || V == n - 2)) return true;

    Int Vk = V;
    for (int i = 1; i < s; i++) {
        if (Vk == 0) return true;
        Vk = mod_n(Vk * Vk - 2, n); // Q = 1
    }
    return false;
}

bool is_almost_extra_strong_lucas_prp(Int n, Int D, Int P, Int Q) {
    (void)D; // unused
    Int dd = n + 1;
    int s = 0;
    while (dd % 2 == 0) { dd /= 2; s++; }
    auto [U, V, __] = lucas_sequence(n, P, Q, dd);
    if (V == 2 || V == n - 2) {
        if (U == 0) return true;
    }

    Int Vk = V;
    for (int i = 0; i < s - 1; i++) {
        if (Vk == 0) return true;
        Vk = mod_n(Vk * Vk - 2, n); // Q = 1
    }
    return false;
}

// ---------------------------------------------------------------------------
// Strong probable prime with arbitrary base
// ---------------------------------------------------------------------------

bool strong_probable_prime_base(Int n, Int base) {
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;

    Int d = n - 1;
    int s = 0;
    while (d % 2 == 0) { d /= 2; s++; }

    Int x = modpow(base, d, n);
    if (x == 1 || x == n - 1) return true;

    Int cur = x;
    for (int i = 0; i < s - 1; i++) {
        cur = (cur * cur) % n;
        if (cur == n - 1) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Quick precheck
// ---------------------------------------------------------------------------

bool precheck(Int n, bool& is_prime) {
    if (n < 2) { is_prime = false; return true; }
    if (n == 2) { is_prime = true; return true; }
    if (n % 2 == 0) { is_prime = false; return true; }
    if (n == 3) { is_prime = true; return true; }
    if (n % 3 == 0) { is_prime = false; return true; }
    if (is_square(n)) { is_prime = false; return true; }
    return false; // undecided
}

// ---------------------------------------------------------------------------
// BPSW variants
// ---------------------------------------------------------------------------

bool bpsw_standard(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    auto [D, P, Q] = selfridge_method_a(n);
    if (P == 0) return false;
    return is_strong_lucas_prp(n, D, P, Q);
}

bool bpsw_regular_lucas(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    auto [D, P, Q] = selfridge_method_a(n);
    if (P == 0) return false;
    return is_regular_lucas_prp(n, D, P, Q);
}

bool bpsw_extra_strong(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    auto [D, P, Q] = extra_strong_params(n);
    if (P == 0) return false;
    return is_extra_strong_lucas_prp(n, D, P, Q);
}

bool bpsw_almost_extra_strong(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    auto [D, P, Q] = extra_strong_params(n);
    if (P == 0) return false;
    return is_almost_extra_strong_lucas_prp(n, D, P, Q);
}

bool bfw_enhanced_bpsw(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    auto [D, P, Q] = selfridge_method_a_star(n);
    if (P == 0) return false;
    if (!is_strong_lucas_prp(n, D, P, Q)) return false;
    // Step 4: Lucas-V: V_{n+1} ≡ 2*Q (mod n)
    auto [_, V, __] = lucas_sequence(n, P, Q, n + 1);
    if (V != mod_n(Int(2) * Q, n)) return false;
    // Step 5: Euler criterion: Q^((n+1)/2) ≡ Q * (Q/n) (mod n)
    Int lhs = modpow(Q, (n + 1) / 2, n);
    Int j = jacobi_symbol(Q, n);
    Int rhs = mod_n(Q * Int(j), n);
    if (lhs != rhs) return false;
    return true;
}

bool bpsw21(Int n) { return bfw_enhanced_bpsw(n); }
bool accelerated_bfw(Int n) { return bfw_enhanced_bpsw(n); }

bool superbfpsw(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    auto [D, P, Q] = selfridge_method_a(n);
    if (P == 0) return false;
    Int t = (D < 0 ? -D : D) % n;
    Int base = (t > 1) ? t : 5;
    if (!strong_probable_prime_base(n, base)) return false;
    return is_strong_lucas_prp(n, D, P, Q);
}

bool fibonacci_bpsw(Int n) {
    bool decided;
    if (precheck(n, decided)) return decided;
    if (!is_strong_probable_prime_base2(n)) return false;
    Int n_mod_5 = mod_n(n, 5);
    if (n_mod_5 == 2 || n_mod_5 == 3) {
        // Fibonacci test: F_{n+1} ≡ 0 (mod n), P=1, Q=-1
        auto [U, _, __] = lucas_sequence(n, 1, -1, n + 1);
        return U == 0;
    } else {
        auto [D, P, Q] = selfridge_method_a(n);
        if (P == 0) return false;
        return is_regular_lucas_prp(n, D, P, Q);
    }
}

// ---------------------------------------------------------------------------
// Variant registry
// ---------------------------------------------------------------------------

const std::vector<VariantEntry> VARIANTS = {
    {1, "Regular Lucas BPSW", bpsw_regular_lucas},
    {2, "Standard Strong BPSW", bpsw_standard},
    {3, "Extra Strong Lucas BPSW", bpsw_extra_strong},
    {4, "Almost Extra Strong Lucas BPSW", bpsw_almost_extra_strong},
    {5, "BFW / Enhanced BPSW (2021)", bfw_enhanced_bpsw},
    {6, "BPSW21", bpsw21},
    {7, "Accelerated BFW", accelerated_bfw},
    {8, "SuperBFPSW (Hamburg 2025)", superbfpsw},
    {9, "Fibonacci-BPSW", fibonacci_bpsw},
};
