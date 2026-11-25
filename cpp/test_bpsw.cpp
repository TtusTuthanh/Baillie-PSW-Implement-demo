#include "bpsw.h"
#include <cassert>
#include <cstdio>
#include <vector>

// Convert unsigned long long to Int
inline Int to_int(unsigned long long v) { return (Int)v; }

// Simple trial division primality test (for verifying small numbers)
bool is_prime_trial(Int n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (Int i = 3; i * i <= n; i += 2) {
        if (n % i == 0) return false;
    }
    return true;
}

int main() {
    int passed = 0, failed = 0;

    auto check = [&](bool cond, const char* msg) {
        if (cond) { passed++; }
        else { failed++; printf("FAIL: %s\n", msg); }
    };

    // ---- Small primes ----
    unsigned long long small_primes[] = {
        2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
        73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
        157,163,167,173,179,181,191,193,197,199
    };
    for (auto p : small_primes) {
        Int pp = to_int(p);
        for (const auto& v : VARIANTS) {
            check(v.func(pp), "variant rejected prime");
        }
    }

    // ---- Small composites ----
    unsigned long long small_composites[] = {
        4,6,8,9,10,12,14,15,16,18,20,21,22,24,25,26,27,28,30,32,
        33,34,35,36,38,39,40,42,44,45,46,48,49,50,51,52,54,55,56,57,
        58,60,62,63,64,65,66,68,69,70,72,74,75,76,77,78,80,81,82,84,
        85,86,87,88,90,91,92,93,94,95,96,98,99,100
    };
    for (auto c : small_composites) {
        Int cc = to_int(c);
        for (const auto& v : VARIANTS) {
            check(!v.func(cc), "variant accepted composite");
        }
    }

    // ---- Strong PSP base 2 ----
    unsigned long long psp[] = {2047, 3277, 4033, 4681, 8321, 15841, 29341};
    for (auto n : psp) {
        Int nn = to_int(n);
        for (const auto& v : VARIANTS) {
            check(!v.func(nn), "variant accepted strong PSP base 2");
        }
    }

    // ---- Large primes ----
    Int large_primes[] = {
        to_int(1000003ULL),
        to_int(1000033ULL),
        to_int(1000000007ULL),
        to_int(1000000009ULL),
    };
    for (auto p : large_primes) {
        for (const auto& v : VARIANTS) {
            check(v.func(p), "variant rejected large prime");
        }
    }

    // ---- Large composites ----
    Int large_composites[] = {
        to_int(1000001ULL),
        to_int(1000000001ULL),
    };
    for (auto c : large_composites) {
        for (const auto& v : VARIANTS) {
            check(!v.func(c), "variant accepted large composite");
        }
    }

    // ---- Edge cases ----
    for (const auto& v : VARIANTS) {
        check(!v.func(0), "zero");
        check(!v.func((Int)(-5)), "negative");
        check(v.func(2), "prime 2");
        check(!v.func(4), "composite 4");
        check(!v.func(9), "square 9");
        check(!v.func(25), "square 25");
    }

    // ---- Variant agreement on 2..2000 ----
    for (Int n = 2; n <= 2000; n = n + 1) {
        bool first = VARIANTS[0].func(n);
        for (size_t i = 1; i < VARIANTS.size(); i++) {
            check(VARIANTS[i].func(n) == first, "variants disagree");
        }
        check(first == is_prime_trial(n), "trial division mismatch");
    }

    // ---- Jacobi symbol ----
    check(jacobi_symbol(5, 13) == -1, "jacobi(5,13)");
    check(jacobi_symbol(-7, 13) == -1, "jacobi(-7,13)");
    check(jacobi_symbol(1, 7) == 1, "jacobi(1,7)");
    check(jacobi_symbol(0, 7) == 0, "jacobi(0,7)");
    check(jacobi_symbol(2, 7) == 1, "jacobi(2,7)");
    check(jacobi_symbol(2, 5) == -1, "jacobi(2,5)");

    // ---- Lucas sequence ----
    auto [U10, V10, Q10] = lucas_sequence(1001, 1, -1, 10);
    check(U10 == 55, "lucas F_10");
    check(V10 == 123, "lucas L_10");
    check(Q10 == 1, "lucas Q^10");

    // ---- Selfridge Method A ----
    for (auto p : small_primes) {
        if (p < 5) continue;
        Int pp = to_int(p);
        auto [D, PP, QQ] = selfridge_method_a(pp);
        check(PP == 1, "selfridge P=1");
        check(jacobi_symbol(D, pp) == -1, "selfridge jacobi=-1");
        check(QQ == (1 - D) / 4, "selfridge Q formula");
    }

    // ---- Extra Strong params ----
    for (auto p : small_primes) {
        if (p < 5) continue;
        Int pp = to_int(p);
        auto [D, PP, QQ] = extra_strong_params(pp);
        check(QQ == 1, "extra Q=1");
        check(PP > 2, "extra P>=3");
        check(jacobi_symbol(D, pp) == -1, "extra jacobi=-1");
        check(D == PP * PP - 4, "extra D=P^2-4");
    }

    printf("\nResults: %d passed, %d failed\n", passed, failed);
    return failed > 0 ? 1 : 0;
}
