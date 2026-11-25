"""Baillie-PSW primality test — all documented variants.

Variants implemented:
  1. Regular Lucas BPSW             — MR base 2 + Regular Lucas (U_{n+1} ≡ 0)
  2. Standard Strong BPSW           — MR base 2 + Strong Lucas-Selfridge (Method A)
  3. Extra Strong Lucas BPSW        — MR base 2 + Extra Strong Lucas (Q=1)
  4. Almost Extra Strong Lucas BPSW — MR base 2 + Almost Extra Strong Lucas
  5. BFW / Enhanced BPSW (2021)     — MR + Strong Lucas + Lucas-V + Euler criterion
  6. BPSW21                         — Same Lucas checks as BFW 2021
  7. Accelerated BFW                — Same conditions as BFW, optimized constants
  8. SuperBFPSW (Hamburg 2025)      — Approximation (paper behind Cloudflare)
  9. Fibonacci-BPSW                 — MR + Fibonacci when n ≡ 2,3 mod 5
"""

from .variants import (
    bpsw_standard,
    bpsw_regular_lucas,
    bpsw_extra_strong,
    bpsw_almost_extra_strong,
    is_bfw_enhanced_bpsw,
    is_bpsw21,
    is_accelerated_bfw,
    is_superbfpsw,
    is_fibonacci_bpsw,
    is_regular_lucas_prp,
    is_strong_lucas_prp,
    is_extra_strong_lucas_prp,
    is_almost_extra_strong_lucas_prp,
    VARIANTS,
)
from .core import (
    jacobi_symbol,
    is_square,
    is_strong_probable_prime_base2,
    lucas_sequence,
    selfridge_method_a,
    selfridge_method_a_star,
    extra_strong_params,
)

__all__ = [
    "bpsw_standard",
    "bpsw_regular_lucas",
    "bpsw_extra_strong",
    "bpsw_almost_extra_strong",
    "is_bfw_enhanced_bpsw",
    "is_bpsw21",
    "is_accelerated_bfw",
    "is_superbfpsw",
    "is_fibonacci_bpsw",
    "is_regular_lucas_prp",
    "is_strong_lucas_prp",
    "is_extra_strong_lucas_prp",
    "is_almost_extra_strong_lucas_prp",
    "VARIANTS",
    "jacobi_symbol",
    "is_square",
    "is_strong_probable_prime_base2",
    "lucas_sequence",
    "selfridge_method_a",
    "selfridge_method_a_star",
    "extra_strong_params",
]
