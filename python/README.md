# Python — Baillie-PSW (BPSW) Primality Test

All 9 documented BPSW variants implemented in pure Python (arbitrary-precision `int`).

## Variants

| # | Variant | Lucas test | Parameters |
|---|---------|-----------|------------|
| 1 | Regular Lucas BPSW | Regular (U_{n+1} ≡ 0) | Selfridge Method A |
| 2 | Standard Strong BPSW | Strong Lucas | Selfridge Method A |
| 3 | Extra Strong Lucas BPSW | Extra Strong Lucas | Q=1, find P |
| 4 | Almost Extra Strong Lucas BPSW | Almost Extra Strong | Q=1, find P (Go method C) |
| 5 | BFW / Enhanced BPSW (2021) | Strong Lucas + Lucas-V + Euler | Method A* |
| 6 | BPSW21 | Strong Lucas + Lucas-V + Euler | Method A* |
| 7 | Accelerated BFW | Same as BFW (optimized constants) | Method A* |
| 8 | SuperBFPSW (Hamburg 2025) | Strong Lucas + Strong Fermat base T | Method A (approx.) |
| 9 | Fibonacci-BPSW | Fibonacci when n≡2,3 mod 5 | P=1, Q=−1 |

## Installation

```bash
cd python
uv venv
uv sync
```

## Usage

```python
from bpsw import bpsw_standard, VARIANTS, is_bfw_enhanced_bpsw

# Standard BPSW (recommended)
bpsw_standard(1000000007)  # True

# Run all variants on a number
for num, (name, func) in VARIANTS.items():
    print(f"{name}: {func(1000000007)}")

# BFW 2021 Enhanced
is_bfw_enhanced_bpsw(1000000007)  # True
```

## Testing

```bash
cd python
uv run pytest tests/ -v
```

## Notes

- Variant 8 (SuperBFPSW) is approximate: the original paper (ePrint 2025/2083) was behind Cloudflare and could not be accessed. Only the abstract was obtained.
- All variants use Python's built-in arbitrary-precision integers — no external dependencies required.
- This is a probable-prime test, not a primality proof. For proof, combine with ECPP or APR-CL.
