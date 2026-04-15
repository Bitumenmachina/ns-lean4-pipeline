#!/usr/bin/env python3
"""S1: Generate Lean 4 definitions from params.json.

Converts floating-point experimental data to exact rational arithmetic
for Lean 4 formalization. Uses limit_denominator(10^8) for clean fractions.

Output: generated/NSParams.lean
"""

import json
from fractions import Fraction
from pathlib import Path

PARAMS = Path(__file__).parent.parent / "params.json"
OUT = Path(__file__).parent.parent / "generated" / "NSParams.lean"


def to_rat(x, denom_limit=10**8):
    """Convert float to Lean rational literal string."""
    f = Fraction(x).limit_denominator(denom_limit)
    if f.denominator == 1:
        return f"({f.numerator} : ℚ)"
    return f"({f.numerator} / {f.denominator} : ℚ)"


def main():
    p = json.load(PARAMS.open())
    OUT.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines.append("import Mathlib.Data.Rat.Defs")
    lines.append("import Mathlib.Algebra.Order.Field.Rat")
    lines.append("")
    lines.append("/-!")
    lines.append("  NSParams.lean — Auto-generated from V5 experimental data")
    lines.append(f"  Source: params.json ({PARAMS.stat().st_size} bytes)")
    lines.append(f"  Databases: {', '.join(p['_meta']['source_dbs'])}")
    lines.append("")
    lines.append("  K = Λ · η  where  η = (ν³/ε)^(1/4)  (Kolmogorov scale)")
    lines.append("  All values are exact rationals via limit_denominator(10^8).")
    lines.append("-/")
    lines.append("")
    lines.append("namespace NS")
    lines.append("")

    # ── Solver constants ─────────────────────────────────────────────
    lines.append("-- Solver parameters")
    lines.append(f"noncomputable def nu_Re6400 : ℚ := {to_rat(p['eta']['nu_Re6400'])}")
    lines.append("")

    # ── Key scalar values ────────────────────────────────────────────
    lines.append("-- Key measurements at Re=6400")
    lines.append(f"def Lambda_Re6400_t20 : ℚ := {to_rat(p['Lambda']['Re6400_N128_t20'])}")
    lines.append(f"def Lambda_Re6400_t40 : ℚ := {to_rat(p['Lambda']['Re6400_N128_t40'])}")
    lines.append(f"def Lambda_Re6400_N256_t20 : ℚ := {to_rat(p['Lambda']['Re6400_N256_t20'])}")
    lines.append(f"def eta_Re6400_t40 : ℚ := {to_rat(p['eta']['Re6400_t40'])}")
    lines.append(f"def epsilon_Re6400_t40 : ℚ := {to_rat(p['eta']['dissipation_Re6400_t40'])}")
    lines.append(f"def K_Re6400_t40 : ℚ := {to_rat(p['K']['Re6400_t40'])}")
    lines.append(f"def K_Re6400_t20 : ℚ := {to_rat(p['K']['Re6400_t20'])}")
    lines.append("")

    # ── Enstrophy ────────────────────────────────────────────────────
    lines.append("-- Enstrophy (monotonically decreasing)")
    lines.append(f"def Omega_Re6400_t0 : ℚ := {to_rat(p['enstrophy']['Re6400_t0'])}")
    lines.append(f"def Omega_Re6400_t20 : ℚ := {to_rat(p['enstrophy']['Re6400_t20'])}")
    lines.append(f"def Omega_Re6400_t40 : ℚ := {to_rat(p['enstrophy']['Re6400_t40'])}")
    lines.append("")

    # ── Convergence ratios ───────────────────────────────────────────
    lines.append("-- Resolution convergence ratios (must be < 1.3 and trending → 1)")
    conv = p['convergence']['Re6400_t20']
    lines.append(f"def conv_ratio_128_64 : ℚ := {to_rat(conv['ratio_128_64'])}")
    lines.append(f"def conv_ratio_256_128 : ℚ := {to_rat(conv['ratio_256_128'])}")
    lines.append("")

    # ── Onset table ──────────────────────────────────────────────────
    lines.append("-- Onset table: (Re, Lambda_final, final/trough)")
    lines.append("-- Monotonically increasing final/trough with Re")
    lines.append("def onset_table : List (ℚ × ℚ × ℚ) :=")
    entries = []
    for row in p['onset_table']:
        re = to_rat(row['Re'])
        lf = to_rat(row['Lambda_final'])
        ft = to_rat(row['final_over_trough'])
        entries.append(f"  ({re}, {lf}, {ft})")
    lines.append(",\n".join(entries).replace(entries[0], "  [" + entries[0].strip(), 1) + "]")
    lines.append("")

    # ── Staircase growth ─────────────────────────────────────────────
    lines.append("-- Staircase growth at Re=6400: (t, Lambda, K=Lambda*eta)")
    lines.append("-- Pattern: oscillatory ratchet, each peak higher")
    lines.append("def staircase : List (ℚ × ℚ × ℚ) :=")
    entries = []
    for step in p['staircase_growth']['steps']:
        t = to_rat(step['t'])
        lam = to_rat(step['Lambda'])
        k = to_rat(step['K_Lambda_eta'])
        entries.append(f"  ({t}, {lam}, {k})")
    lines.append(",\n".join(entries).replace(entries[0], "  [" + entries[0].strip(), 1) + "]")
    lines.append("")

    # ── K relationships across Re ────────────────────────────────────
    lines.append("-- K = Λ·η at t=20 for each Re (K decreases with Re at fixed t)")
    lines.append("def K_vs_Re : List (ℚ × ℚ) :=")
    entries = []
    for row in p['K']['relationships']:
        re = to_rat(row['Re'])
        k = to_rat(row['K'])
        entries.append(f"  ({re}, {k})")
    lines.append(",\n".join(entries).replace(entries[0], "  [" + entries[0].strip(), 1) + "]")
    lines.append("")

    # ── Re bifurcation bounds ────────────────────────────────────────
    lines.append("-- Bifurcation regime boundaries (continuous, not sharp)")
    bif = p['Re_bifurcation']
    lines.append(f"def Re_regime1_upper : ℚ := {to_rat(bif['regime_1_upper'])}")
    lines.append(f"def Re_regime2_lower : ℚ := {to_rat(bif['regime_2_range'][0])}")
    lines.append(f"def Re_regime2_upper : ℚ := {to_rat(bif['regime_2_range'][1])}")
    lines.append(f"def Re_regime3_lower : ℚ := {to_rat(bif['regime_3_lower'])}")
    lines.append("")

    # ── Kill metric definition ───────────────────────────────────────
    lines.append("/-- The kill metric K = Λ · η. Dimensionless.")
    lines.append("    Direction field gradient (Λ) scaled by Kolmogorov length (η).")
    lines.append("    K > 1 means direction varies faster than dissipation scale. -/")
    lines.append("noncomputable def K_metric (Lambda eta : ℚ) : ℚ := Lambda * eta")
    lines.append("")

    lines.append("end NS")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT}  ({OUT.stat().st_size} bytes)")
    n_defs = sum(1 for l in lines if l.strip().startswith("def ") or l.strip().startswith("noncomputable def "))
    print(f"  {n_defs} definitions")


if __name__ == "__main__":
    main()
