#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

// Use __int128 for arbitrary-precision arithmetic (GCC/Clang extension)
// Supports integers up to ~2^127 (~1.7e38), sufficient for most practical tests.
using Int = __int128;
using UInt = unsigned __int128;

// ---------------------------------------------------------------------------
// Core mathematical primitives
// ---------------------------------------------------------------------------

// Jacobi symbol (a/n). n must be odd and positive.
int jacobi_symbol(Int a, Int n);

// Check if n is a perfect square.
bool is_square(Int n);

// Strong probable prime test base 2 (Miller-Rabin).
bool is_strong_probable_prime_base2(Int n);

// Compute (U_k mod n, V_k mod n, Q^k mod n) using binary ladder.
std::tuple<Int, Int, Int> lucas_sequence(Int n, Int P, Int Q, Int k);

// Selfridge Method A: scan D in 5, -7, 9, -11, ... until Jacobi(D,n) = -1.
// Returns (D, P=1, Q=(1-D)/4). Returns (D, 0, 0) if factor found.
std::tuple<Int, Int, Int> selfridge_method_a(Int n);

// Selfridge Method A*: same as A, but if Q == -1 (D == 5), set P = Q = 5.
std::tuple<Int, Int, Int> selfridge_method_a_star(Int n);

// Extra Strong parameter selection: Q=1, find P from 3 until Jacobi(P^2-4,n) = -1.
std::tuple<Int, Int, Int> extra_strong_params(Int n);

// ---------------------------------------------------------------------------
// Lucas probable prime tests
// ---------------------------------------------------------------------------

bool is_regular_lucas_prp(Int n, Int D, Int P, Int Q);
bool is_strong_lucas_prp(Int n, Int D, Int P, Int Q);
bool is_extra_strong_lucas_prp(Int n, Int D, Int P, Int Q);
bool is_almost_extra_strong_lucas_prp(Int n, Int D, Int P, Int Q);

// ---------------------------------------------------------------------------
// BPSW variants (all return true = probable prime, false = composite)
// ---------------------------------------------------------------------------

// Variant 2: Standard Strong BPSW — MR base 2 + Strong Lucas-Selfridge (Method A).
bool bpsw_standard(Int n);

// Variant 1: Regular Lucas BPSW — MR base 2 + Regular Lucas.
bool bpsw_regular_lucas(Int n);

// Variant 3: Extra Strong Lucas BPSW — MR base 2 + Extra Strong Lucas.
bool bpsw_extra_strong(Int n);

// Variant 4: Almost Extra Strong Lucas BPSW — MR base 2 + Almost Extra Strong.
bool bpsw_almost_extra_strong(Int n);

// Variant 5: BFW / Enhanced BPSW (2021).
bool bfw_enhanced_bpsw(Int n);

// Variant 6: BPSW21 — same Lucas checks as BFW 2021.
bool bpsw21(Int n);

// Variant 7: Accelerated BFW — same conditions as BFW.
bool accelerated_bfw(Int n);

// Variant 8: SuperBFPSW (Hamburg 2025) — approximation (paper behind Cloudflare).
bool superbfpsw(Int n);

// Variant 9: Fibonacci-BPSW — MR base 2 + Fibonacci when n ≡ 2,3 mod 5.
bool fibonacci_bpsw(Int n);

// ---------------------------------------------------------------------------
// Variant registry
// ---------------------------------------------------------------------------

struct VariantEntry {
    int id;
    const char* name;
    bool (*func)(Int);
};

extern const std::vector<VariantEntry> VARIANTS;
