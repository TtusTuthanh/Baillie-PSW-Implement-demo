# Status: BPSW Implementation — COMPLETE

## Completed — All Three Languages

### Python (uv) — `python/`
- core.py, variants.py, __init__.py, tests/test_bpsw.py, pyproject.toml, README.md
- **352 tests PASSED** (pytest)
- Tested: small primes, composites, strong PSP base 2, large primes up to 10^18, variant agreement 2..2000 vs trial division

### Rust — `rust/`
- Cargo.toml, src/lib.rs (with unit tests), uses num-bigint 0.4
- **12 tests PASSED** (cargo test)
- Same test coverage as Python

### C++ — `cpp/`
- CMakeLists.txt, bpsw.h, bpsw.cpp, test_bpsw.cpp, uses __int128 (GCC/Clang extension)
- **19,559 assertions PASSED** (ctest)
- Same test coverage plus larger range (2..2000)

## 9 Variants Implemented (all languages)
| # | Variant | Lucas Test | Parameters |
|---|---------|-----------|------------|
| 1 | Regular Lucas BPSW | Regular (U_{n+1}≡0) | Selfridge Method A |
| 2 | Standard Strong BPSW | Strong Lucas | Selfridge Method A |
| 3 | Extra Strong Lucas BPSW | Extra Strong Lucas | Q=1, find P |
| 4 | Almost Extra Strong Lucas BPSW | Almost Extra Strong | Q=1, find P (Go method C) |
| 5 | BFW / Enhanced BPSW (2021) | Strong Lucas + Lucas-V + Euler | Method A* |
| 6 | BPSW21 | Same as BFW 2021 | Method A* |
| 7 | Accelerated BFW | Same conditions as BFW | Method A* |
| 8 | SuperBFPSW (Hamburg 2025) | Strong Lucas + Strong Fermat base T | Method A (approx.) |
| 9 | Fibonacci-BPSW | Fibonacci when n≡2,3 mod 5 | P=1, Q=−1 |

## Gap
- Variant 8 (SuperBFPSW): original paper (ePrint 2025/2083) behind Cloudflare. Implemented as approximation based on abstract, clearly marked in all three language implementations.

## Reference Material (saved in `docs/`)
- sympy-primetest.py, go-prime.decoded.go, bfw-extracted.md (BFW 2021 from arXiv HTML)
- chen-greene-2003 (reconstructed), hamburg-superbfpsw-2025.pdf (Cloudflare block)

## Cross-Language Verification
All three implementations agree on:
- Primes 2..199 accepted by all variants
- Composites up to 100 rejected by all variants
- Strong PSP base 2 (2047..29341) rejected by all variants
- Large primes 10^6..10^18 accepted
- Variant agreement on 2..2000 matches trial division exactly
