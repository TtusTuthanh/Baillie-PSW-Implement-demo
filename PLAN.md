# Implementation Plan: Baillie-PSW (BPSW) — All Variants

## Overview
Implement all documented BPSW variants in three languages (Python with uv, Rust, C++), with shared mathematical primitives (big-int modular arithmetic, Jacobi symbol, Miller-Rabin base 2, Lucas sequences U/V, parameter selection methods A–D and R1–R2), plus tests and cross-language verification.

## Variants to Implement (9 total)
1. Regular Lucas BPSW — MR base 2 + Regular Lucas (U_{n+1} ≡ 0)
2. Standard Strong BPSW — MR base 2 + Strong Lucas-Selfridge (Method A)
3. Extra Strong Lucas BPSW — MR base 2 + Extra Strong Lucas (Q=1, find P)
4. Almost Extra Strong Lucas BPSW — MR base 2 + Almost Extra Strong (skip U_d=0)
5. BFW / Enhanced BPSW — MR base 2 + Strong Lucas + Lucas-V + Euler criterion
6. BPSW21 — Strong Lucas + Lucas-V + Euler criterion (LucasCheck mode)
7. Accelerated BFW — Optimized BFW implementation
8. SuperBFPSW — Hamburg's improved variant (Strong Fermat base 2 + base T + accelerated Lucas)
9. Fibonacci-BPSW — MR base 2 + Fibonacci test when n ≡ 2,3 mod 5

## Architecture Decisions
- Python: built-in `int` (arbitrary precision), `uv` for project mgmt
- Rust: `num-bigint` + `num-integer` crates
- C++: `boost::multiprecision::cpp_int` (or `__int128` fallback for < 2^128)
- All variants share a common `core` module; variant-specific logic in separate modules
- Parameter selection: Selfridge Method A as default (D in 5,-7,9,-11,..., P=1, Q=(1-D)/4)

## Task List

### Phase 1: Foundation (config + core math)
- [ ] Task 1: Python project scaffold (pyproject.toml, uv, package layout)
- [ ] Task 2: Rust project scaffold (Cargo.toml, dependencies)
- [ ] Task 3: C++ project scaffold (CMakeLists.txt, build config)
- [ ] Task 4: Python core — jacobi_symbol, modpow, miller_rabin_base2, lucas_uv, selfridge_method_A
- [ ] Task 5: Rust core — same primitives
- [ ] Task 6: C++ core — same primitives

### Checkpoint: Foundation
- [ ] All three projects build/compile
- [ ] Core math functions tested against known values

### Phase 2: Lucas sequence variants
- [ ] Task 7: Python — regular_lucas_prp, strong_lucas_prp, extra_strong_lucas_prp, almost_extra_strong_lucas_prp
- [ ] Task 8: Rust — same Lucas variants
- [ ] Task 9: C++ — same Lucas variants

### Checkpoint: Lucas Variants
- [ ] Each Lucas PRP tested against known primes and pseudoprimes

### Phase 3: BPSW algorithm variants
- [ ] Task 10: Python — all 9 BPSW variants
- [ ] Task 11: Rust — all 9 BPSW variants
- [ ] Task 12: C++ — all 9 BPSW variants

### Checkpoint: BPSW Variants
- [ ] Each variant returns correct probable-prime/composite for test vectors

### Phase 4: Tests & verification
- [ ] Task 13: Python tests (pytest) — known primes, composites, BPSW pseudoprimes
- [ ] Task 14: Rust tests — same test vectors
- [ ] Task 15: C++ tests — same test vectors
- [ ] Task 16: Cross-language verification — all three agree on test vectors
- [ ] Task 17: Benchmark comparison (optional)

### Checkpoint: Complete
- [ ] All tests pass in all three languages
- [ ] Results match across languages
- [ ] README with usage and variant comparison table

## Risks and Mitigations
| Risk | Impact | Mitigation |
|------|--------|------------|
| Big-int performance in C++ without GMP | High | Use boost::multiprecision or __int128 for <2^128 |
| Subtle bugs in Lucas sequence implementation | High | Cross-verify against SymPy/Go reference implementations |
| Parameter selection edge cases (Q=-1, Q=1, square D) | Medium | Handle all edge cases from BFW appendix |
| No git repo → no worktrees → one writer at a time | Medium | Sequence language implementations |

## Open Questions
- None — spec is complete from user document
