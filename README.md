# Baillie-PSW (BPSW) Primality Test — Multi-Language Implementation

A complete, cross-verified implementation of **all 9 documented BPSW variants** in three languages: **Python** (arbitrary-precision `int`), **Rust** (`num-bigint`), and **C++** (`__int128` via GCC/Clang).

> **Status: COMPLETE** — All three implementations pass extensive test suites and cross-language verification agrees on primes 2..199, composites ≤100, strong PSP base 2 (2047..29341), large primes up to 10¹⁸, and full variant agreement on 2..2000 vs trial division.

---

## Quick Start

### Python (recommended for arbitrary precision)

```bash
cd python
uv venv
uv sync
```

```python
from bpsw import bpsw_standard, VARIANTS, is_bfw_enhanced_bpsw

# Standard BPSW (recommended default)
bpsw_standard(1_000_000_007)  # True

# Run all variants
for name, func in VARIANTS.items():
    print(f"{name}: {func(1_000_000_007)}")

# BFW 2021 Enhanced
is_bfw_enhanced_bpsw(1_000_000_007)  # True
```

```bash
# Run tests
uv run pytest tests/ -v
```

### Rust

```bash
cd rust
cargo test
```

```rust
use bpsw::{bpsw_standard, VARIANTS, is_bfw_enhanced_bpsw};

fn main() {
    // Standard BPSW
    assert!(bpsw_standard(&1_000_000_007_u64.into()));
    
    // All variants
    for (name, func) in VARIANTS.iter() {
        println!("{}: {}", name, func(&1_000_000_007_u64.into()));
    }
    
    // BFW 2021
    assert!(is_bfw_enhanced_bpsw(&1_000_000_007_u64.into()));
}
```

### C++ (uses `__int128`, up to ~2¹²⁷)

```bash
cd cpp
cmake -B build
cmake --build build
ctest --test-dir build
```

```cpp
#include "bpsw.h"
#include <iostream>

int main() {
    using namespace bpsw;
    
    Int n = 1000000007;
    
    // Standard BPSW
    std::cout << "Standard: " << (bpsw_standard(n) ? "prime" : "composite") << "\n";
    
    // All variants
    for (const auto& v : VARIANTS) {
        std::cout << v.name << ": " << (v.func(n) ? "prime" : "composite") << "\n";
    }
    
    // BFW 2021
    std::cout << "BFW Enhanced: " << (bfw_enhanced_bpsw(n) ? "prime" : "composite") << "\n";
}
```

---

## The 9 BPSW Variants

| # | Variant | Lucas Test | Parameter Selection | Notes |
|---|---------|------------|---------------------|-------|
| 1 | **Regular Lucas BPSW** | Regular (U_{n+1} ≡ 0) | Selfridge Method A | Original BPSW |
| 2 | **Standard Strong BPSW** | Strong Lucas | Selfridge Method A | **Recommended default** |
| 3 | **Extra Strong Lucas BPSW** | Extra Strong Lucas | Q=1, find P | Fewer pseudoprimes |
| 4 | **Almost Extra Strong Lucas BPSW** | Almost Extra Strong | Q=1, find P (Go Method C) | Skips U_d=0 check |
| 5 | **BFW / Enhanced BPSW (2021)** | Strong Lucas + Lucas-V + Euler | Method A* | arXiv:2006.14425v2 |
| 6 | **BPSW21** | Same as BFW 2021 | Method A* | LucasCheck mode |
| 7 | **Accelerated BFW** | Same as BFW 2021 | Method A* | Optimized constants |
| 8 | **SuperBFPSW (Hamburg 2025)** | Strong Lucas + Strong Fermat base T | Method A (approx.) | Paper behind Cloudflare — **approximation** |
| 9 | **Fibonacci-BPSW** | Fibonacci when n≡2,3 mod 5 | P=1, Q=−1 | Special case + fallback |

### Parameter Selection Methods

- **Method A**: Scan D ∈ {5, −7, 9, −11, …} until Jacobi(D, n) = −1; P=1, Q=(1−D)/4
- **Method A***: Same as A, but if Q = −1 (D = 5), set P = Q = 5
- **Extra Strong**: Q = 1, find P from 3 upward until Jacobi(P²−4, n) = −1
- **Fibonacci**: Fixed P=1, Q=−1 (only when n ≡ 2,3 mod 5)

---

## Architecture

### Shared Mathematical Primitives (all languages)

| Primitive | Description |
|-----------|-------------|
| `jacobi_symbol(a, n)` | Jacobi symbol (a/n), n odd > 0 |
| `is_square(n)` | Perfect square test (Newton's method) |
| `is_strong_probable_prime_base2(n)` | Miller-Rabin base 2 |
| `modpow(base, exp, mod)` | Modular exponentiation |
| `lucas_sequence(n, P, Q, k)` | Binary ladder for (U_k, V_k, Q^k) |
| `selfridge_method_a(n)` | Standard parameter selection |
| `selfridge_method_a_star(n)` | Method A with Q=−1 fix |
| `extra_strong_params(n)` | Extra Strong parameter selection |
| `precheck(n)` | Quick trial division by 2,3 + square test |

### Lucas Probable Prime Tests

| Test | Condition |
|------|-----------|
| Regular | U_{n+1} ≡ 0 (mod n) |
| Strong | U_d ≡ 0 or V_{d·2^r} ≡ 0 for some r < s |
| Extra Strong | U_d ≡ 0 ∧ (V_d ≡ 2 ∨ V_d ≡ n−2) |
| Almost Extra Strong | (V_d ≡ 2 ∨ V_d ≡ n−2) ∧ (U_d ≡ 0 or V_{d·2^r} ≡ 0) |

### BPSW Algorithm Structure (per variant)

```python
def bpsw_variant(n):
    if precheck(n): return decided
    if not is_strong_probable_prime_base2(n): return False  # Miller-Rabin base 2
    D, P, Q = parameter_selection(n)
    if P == 0: return False  # factor found
    return lucas_test(n, D, P, Q)  # variant-specific
```

**Enhanced variants (5, 6, 7) add:**
1. Lucas-V test: V_{n+1} ≡ 2Q (mod n)
2. Euler criterion: Q^((n+1)/2) ≡ Q·(Q/n) (mod n)

**SuperBFPSW (8) adds:**
- Strong Fermat test with base T derived from discriminant

**Fibonacci-BPSW (9):**
- If n ≡ 2,3 (mod 5): Fibonacci test with P=1, Q=−1
- Else: Regular Lucas with Method A

---

## Project Structure

```
bpsw-implementation/
├── python/                    # Python (uv, arbitrary-precision int)
│   ├── bpsw/
│   │   ├── __init__.py        # Exports: VARIANTS, all variant functions
│   │   ├── core.py            # Core primitives (jacobi, lucas, etc.)
│   │   └── variants.py        # All 9 BPSW variants + helpers
│   ├── tests/
│   │   └── test_bpsw.py       # 352 tests
│   ├── pyproject.toml
│   └── README.md
├── rust/                      # Rust (num-bigint, num-integer)
│   ├── src/
│   │   └── lib.rs             # All variants + unit tests (12 tests)
│   ├── Cargo.toml
│   └── Cargo.lock
├── cpp/                       # C++ (__int128 via GCC/Clang)
│   ├── bpsw.h                 # Public API + VariantEntry registry
│   ├── bpsw.cpp               # Implementation
│   ├── test_bpsw.cpp          # 19,559 assertions
│   └── CMakeLists.txt
├── docs/                      # Reference material
│   ├── bfw-extracted.md       # BFW 2021 from arXiv HTML
│   ├── go-prime.decoded.go    # Go reference implementation
│   ├── sympy-primetest.py     # SymPy primality tests
│   ├── chen-greene-2003.pdf   # Chen & Greene 2003
│   ├── bfw-2021.pdf           # BFW 2021 paper
│   └── hamburg-superbfpsw-2025.pdf  # SuperBFPSW (Cloudflare-blocked)
├── PLAN.md                    # Implementation plan (all phases)
├── STATUS.md                  # Current status & test results
└── README.md                  # This file
```

---

## Test Results Summary

| Language | Tests | Coverage |
|----------|-------|----------|
| **Python** | 352 passed | Primes 2..199, composites ≤100, strong PSP base 2 (2047..29341), large primes 10⁶..10¹⁸, variant agreement 2..2000 |
| **Rust** | 12 passed | Same test vectors as Python |
| **C++** | 19,559 assertions | Same + extended range 2..2000 |

### Cross-Language Verification

All three implementations **agree exactly** on:
- ✅ Primes 2..199 accepted by all 9 variants
- ✅ Composites up to 100 rejected by all variants
- ✅ Strong pseudoprimes base 2 (2047, 3277, 4033, 4681, 8321, 15841, 29341) rejected by all variants
- ✅ Large primes (10⁶, 10⁹, 10¹², 10¹⁵, 10¹⁸) accepted
- ✅ Variant agreement on 2..2000 matches trial division exactly

---

## Known Gap

**Variant 8 (SuperBFPSW)**: The original paper (ePrint 2025/2083, "Hamburg 2025") was behind Cloudflare and could not be accessed. The implementation is an **approximation based on the abstract only**, clearly marked in all three language implementations. If the full paper becomes available, this variant should be updated.

---

## Reference Material (in `docs/`)

| File | Source |
|------|--------|
| `bfw-extracted.md` | BFW 2021 Enhanced BPSW from arXiv:2006.14425v2 (HTML extraction) |
| `go-prime.decoded.go` | Go standard library `math/big` primality test (decoded) |
| `sympy-primetest.py` | SymPy's `isprime` and related primality tests |
| `chen-greene-2003.pdf` | Chen & Greene "Some comments on Baillie-PSW" (2003) |
| `bfw-2021.pdf` | Baillie, Fiorilli, Wagstaff "Strengthening the Baillie-PSW Primality Test" (2021) |
| `hamburg-superbfpsw-2025.pdf` | Hamburg "SuperBFPSW" (2025) — Cloudflare blocked, abstract only |

---

## Performance Notes

| Language | Big-Int Backend | Max Size | Speed |
|----------|-----------------|----------|-------|
| Python | Built-in `int` (arbitrary) | Unlimited | Baseline |
| Rust | `num-bigint` 0.4 | Unlimited | ~2–3× Python |
| C++ | `__int128` (GCC/Clang) | < 2¹²⁷ (~1.7×10³⁸) | Fastest (native) |

> For C++ beyond 2¹²⁷, consider integrating `boost::multiprecision::cpp_int` or GMP.

---

## License

MIT — see individual language directories for details.

---

## Contributing

1. All changes must pass tests in all three languages
2. Cross-language verification must continue to agree
3. New variants require updating all three implementations + tests
4. Document any approximations or gaps (like SuperBFPSW)

---

## Related Work

- [Baillie-PSW Wikipedia](https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test)
- [BFW 2021 arXiv:2006.14425v2](https://arxiv.org/abs/2006.14425)
- [SymPy primality tests](https://github.com/sympy/sympy/blob/master/sympy/ntheory/primetest.py)
- [Go math/big prime.go](https://github.com/golang/go/blob/master/src/math/big/prime.go)