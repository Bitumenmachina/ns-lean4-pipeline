# Axiom Audit Trail

Each `axiom` in `NSCerts.lean` traces to a specific assertion block in
`z3/all_certs.smt2`. Z3 proved all 10 certificates SAT in QF_LRA.

## Axiom → SMT2 Mapping

| Lean Axiom | Z3 Certificate | SMT2 Block | What It Verifies |
|-----------|----------------|------------|-----------------|
| `K_Re6400_t40_gt_one` | `K_t40_gt_1` | lines 7-11 | K(Re=6400,t=40) = 37.72 > 1 |
| `K_early_lt_one` | `K_crossover` | lines 13-20 | K(t≈5) = 0.87 < 1 |
| `K_late_gt_one` | `K_crossover` | lines 13-20 | K(t≈14) = 1.74 > 1 |
| `conv_ratio_256_128_lt_13` | `convergence_ratio_bounded` | lines 22-28 | ratio(256/128) = 1.0995 < 1.3 |
| `conv_ratios_decreasing` | `convergence_ratios_decreasing` | lines 30-37 | ratio(128/64) > ratio(256/128) |
| `onset_ft_400_lt_800` | `onset_ft_strictly_increasing` | lines 39-53 | f/t monotone: 0.9 < 1.1 |
| `onset_ft_800_lt_1600` | `onset_ft_strictly_increasing` | lines 39-53 | f/t monotone: 1.1 < 1.9 |
| `onset_ft_4800_lt_5600` | `onset_ft_strictly_increasing` | lines 39-53 | f/t monotone: 7.8 < 8.8 |
| `onset_ft_5600_lt_6400` | `onset_ft_strictly_increasing` | lines 39-53 | f/t monotone: 8.8 < 9.8 |
| `peak_K_20_lt_30` | `staircase_peaks_increasing` | lines 55-67 | K peaks: 6.9 < 32.21 |
| `peak_K_30_lt_39` | `staircase_peaks_increasing` | lines 55-67 | K peaks: 32.21 < 36.79 |
| `peak_K_39_lt_40` | `staircase_peaks_increasing` | lines 55-67 | K peaks: 36.79 < 37.72 |
| `Lambda_grows_while_Omega_decays` | `Lambda_grows_Omega_decays` | lines 69-78 | Λ: 392→1950, Ω: 9.99→6.88 |
| `ft_Re400_le` | `subcritical_tame` | lines 80-87 | Re≤800: f/t ≤ 1.1 |
| `ft_Re800_le` | `subcritical_tame` | lines 80-87 | Re≤800: f/t ≤ 1.1 |
| `K_product_close` | `K_def_consistent` | lines 89-97 | Λ·η ≈ K (within 0.5) |
| `regime_boundaries_ordered` | `regime_boundaries_ordered` | lines 99-110 | 0 < 800 < 1600 < 4800 < 5600 |

## Trust Boundaries

### Formally Verified (Lean 4 kernel)
All 13 theorems in NSTheorems.lean. The Lean kernel is the trusted
computing base. No `sorry`, no `native_decide`, no `Decidable` escape hatches.

### Mechanically Verified (Z3 QF_LRA)
The 10 certificates above. QF_LRA is a decidable fragment — Z3's result
is sound. Re-audit: `python3 z3/run_z3.py`.

### Empirical (Not Formally Verifiable)
Two hypotheses in the main theorem (OBL-10):
1. **H1:** Λ(t) → ∞ as t → T (direction crease forms in finite time)
2. **H2:** |ω| > c > 0 at the crease (vorticity doesn't vanish there)

**Status after M4 analysis (2026-04-12):**

H1 is supported: sup|∇ξ|_F scales as N^{0.86 ± 0.05} (R² > 0.96), consistent
with O(1/dx). The raw Λ grows with resolution.

H2 is NOT independently supported at the peak-gradient voxel: |ω| at the
location of sup|∇ξ| drops to <0.001 at all resolutions and decreases with N.
Masked-Λ analysis shows the O(1/dx) scaling collapses at θ ≥ 0.05.

**However:** The product sup(|ω|·|∇ξ|_F) is resolution-independent
(γ = 0.04 ≈ 0, CV < 3%). At the spatial location of this product's maximum
(the interface ring between tubes), both |ω| and |∇ξ| are moderate and scale
inversely. H1 and H2 are coupled through this product invariant — they are
not independently realizable in the current data.

The theorem remains valid as a conditional statement. The data characterizes
the scaling balance at the crease rather than instantiating the antecedents
independently. See `data/crease_invariants/crease_invariants_report.md` and
`OPEN_QUESTIONS.md` for the full analysis and unresolved questions.

Finite computation cannot prove infinite-time divergence.
The theorem is conditional: IF H1 ∧ H2 THEN regularity fails.
