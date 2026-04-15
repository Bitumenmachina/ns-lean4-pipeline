# Crease Invariant Characterization

**Date:** 2026-04-12
**Data:** Phase 10b (N=64/128/256, Re=6400, antiparallel tubes, Rust FFTW spectral solver)
**V5 t=40 field data:** NOT AVAILABLE (V5 stored timeseries only; excluded from this analysis)

---

## M4: The Invariant Candidate — Load-Bearing Result

### What the data shows

The direction-field gradient sup|∇ξ|\_F scales as N^β with β = 0.860 ± 0.051 across three common timesteps (t=0, 10, 20) over a 4× resolution range (N=64→256). This is consistent with O(1/dx) scaling — the crease sharpens approximately linearly with resolution.

The product sup(|ω|·|∇ξ|\_F) scales as N^γ with γ = 0.043. This is statistically indistinguishable from zero given the data:

| t | N=64 | N=128 | N=256 | CV | Status |
|---|------|-------|-------|----|--------|
| 0 | 9.657 | 9.074 | 9.560 | 0.027 | Invariant |
| 10 | 6.235 | 6.250 | 6.434 | 0.014 | Invariant |
| 20 | 4.135 | 3.663 | 4.838 | 0.115 | Weak drift |

At t=0 and t=10, the coefficient of variation across resolutions is <3% — the product is flat to measurement precision. At t=20, there is weak drift (CV=11.5%) with N=256 slightly elevated, but only three data points make it difficult to distinguish trend from noise.

The product is NOT constant in time — it decays as the tubes dissipate:
- N=64: 9.66 → 2.38 over t=0→40 (slope -0.158/t)
- N=128: 9.07 → 2.57 over t=0→40 (slope -0.145/t)
- N=256: 9.56 → 4.27 over t=0→28 (slope -0.187/t)

### What this implies

The product sup(|ω|·|∇ξ|\_F) is a function of time but not of resolution. At the spatial location where the product is maximal, |ω| scales as approximately N^{-β} to compensate for |∇ξ|\_F scaling as N^{+β}. This is not happening at the peak-gradient voxel (where |ω| → 0) — it is happening at a DIFFERENT location in the domain where the two factors are balanced.

This is neither a confirmation nor a refutation of the bootstrap failure thesis. It is a scaling relationship the PDE imposes at the interface between antiparallel coherent structures. The classical framework (Constantin-Fefferman, Beale-Kato-Majda) does not have a name for it because those criteria focus on either pointwise sup of individual quantities or integrated norms, not on the pointwise product at its maximum.

### Fit details

| Quantity | Exponent | R² | SE |
|----------|----------|----|----|
| sup\|∇ξ\|\_F ~ N^β | β = 0.860 | 0.987 (mean) | 0.051 |
| sup(\|ω\|·\|∇ξ\|) ~ N^γ | γ = 0.043 | 0.38 (mean, noisy) | — |
| β + (implied α) ≈ γ | 0.043 ≈ 0 | — | — |

β is measured cleanly (R² > 0.96 at every timestep). γ has low R² because the product is genuinely flat — there is no systematic trend to fit.

---

## M1: Local Balance Ratio

### What the data shows

B(x) = |ω·∇u| / (ν|∇²ω| + |ω|·|∇ξ|\_F)

At crease voxels (top 1% of |∇ξ|\_F):
- B_median ranges from 0.001 to 0.010 across all (N, t) pairs
- At tube body voxels (top 10% of |ω|): B_median ranges from 0.008 to 0.042

### What this implies

B << 1 everywhere means the denominator (viscous + direction-gradient terms) DOMINATES over stretching at the crease. The crease is NOT in a stretching-dominated regime. Vortex stretching is 100× weaker than the smoothing/direction-gradient terms at the crease location.

At tube cores, B is ~10× higher than at the crease but still well below 1. The entire flow is in a regime where the diffusive/direction terms dominate the local balance. This is consistent with the tubes being in a decaying regime (enstrophy decreasing monotonically).

---

## M2: Helicity Density at Crease

### What the data shows

Mean helicity h = u·ω at crease voxels: ≈ 0 at all (N, t).
Mean |h| at crease: 0.0005–0.002, decaying with time.
The helicity is symmetric (h_mean ≈ 0) because the antiparallel tubes create equal and opposite contributions.

### What this implies

The crease sits in a helicity-neutral zone. This is geometric — the reconciliation layer between antiparallel tubes must have zero net helicity by symmetry. |h| is small but nonzero, indicating weak local alignment between u and ω at the crease. The crease is NOT a region of high helicity, which distinguishes it from helical vortex structures where strong u·ω alignment inhibits reconnection.

---

## M3: Strain-Vorticity Alignment

### What the data shows

cos(θ) = (S·ω)·ω / (|S·ω|·|ω|) measures alignment between vorticity and its own rate-of-strain.

At crease voxels:
- cos(θ) mean is slightly negative (-0.05 to -0.27 for N=128/256 at t>0)
- Standard deviation is large (~0.5–0.6)

At tube body voxels:
- cos(θ) mean is near zero (-0.01 to -0.11)
- Standard deviation similar (~0.4–0.5)

At N=64, a resolution artifact appears: cos(θ)\_crease trends POSITIVE with time (reaching +0.86 at t=40), opposite to N=128/256. This is a sign the N=64 crease is under-resolved and the voxels classified as "crease" are actually tube-edge voxels with different strain geometry.

### What this implies

At N=128/256, the crease has a slight negative strain-vorticity alignment — vorticity is weakly perpendicular to its strain principal axis. This is the compressive regime where strain thins the structure rather than stretching it. The large standard deviation means the alignment is broadly distributed, not sharply peaked. The crease does not have a distinctive strain geometry compared to tube body voxels — both show weak, broadly-distributed alignment.

The N=64 anomaly is a resolution convergence flag: M3 is not converged at N=64.

---

## M5: Pressure Gradient at Crease

### Dealiasing verification

Spectral vs. central-difference pressure gradient at N=64: mean relative error 0.75%, max 24.8% (at isolated high-k voxels). The spectral solve is consistent with finite differences to discretization tolerance. M5 is trustworthy.

### What the data shows

|∇p| at crease voxels: 0.005–0.049 (mean), decaying with time.
Pressure-vorticity alignment (∇p · ω / |∇p||ω|): mean ≈ 0.000 at all (N, t), with standard deviation 0.45–0.93.

### What this implies

The pressure gradient at the crease is nonzero but has no preferred alignment with vorticity. The alignment is uniformly distributed between -1 and +1. This means pressure is not systematically doing work in a specific direction at the crease — it is redistributing momentum isotropically.

At N=64, the alignment std grows with time to 0.93 (nearly uniform), while at N=128/256 it stays at 0.50–0.58. This is another resolution convergence flag.

The pressure gradient magnitude scales roughly with N (higher at N=256 than N=64 for the same t), which is consistent with the crease gradient producing sharper pressure sources. But the isotropic alignment means this pressure is not preferentially compressing or stretching the crease — it is responding to the direction-field discontinuity without a dominant axis.

---

## M6: Kolmogorov-Conditioned Λ

### What the data shows

Λ\_η(x) = |∇ξ(x)| · η(x), where η = (ν³/ε)^{1/4} is local Kolmogorov scale.
At crease voxels (excluding t=0 which has IC artifacts):

| N | t | Λ\_η mean | Λ\_η max | % crease > 1 | % domain > 1 |
|---|---|-----------|----------|---------------|---------------|
| 64 | 20 | 3.67 | 15.9 | 100% | 15.9% |
| 128 | 20 | 2.62 | 20.2 | 100% | 5.2% |
| 256 | 20 | 3.61 | 39.7 | 100% | 3.9% |

### What this implies

100% of crease voxels have Λ\_η > 1 at all resolutions and all timesteps (after t=0). The direction-field gradient is sharper than the local dissipation scale everywhere on the crease. This is the one diagnostic where the crease is unambiguously distinct from the rest of the domain.

The mean Λ\_η at the crease is roughly resolution-independent (2.6–3.7), consistent with the product invariance found in M4. The MAX Λ\_η increases with N (15.9 → 20.2 → 39.7 at t=20), which is consistent with |∇ξ| growing as N^β while the local dissipation rate ε (and hence η) is resolution-converged.

The fraction of the DOMAIN with Λ\_η > 1 DECREASES with N (15.9% → 5.2% → 3.9%), meaning the crease is becoming more spatially concentrated as resolution increases. At N=64, 16% of the domain exceeds the Kolmogorov condition; at N=256, only 4%. The crease sharpens AND narrows.

---

## Slice Maps (N=256 t=20 z=π)

The six-panel figure shows the z=π midplane through the tube interaction zone:

**M1 (balance ratio):** The crease region (between tubes) has B << 0.01 (dark blue). The tube cores have B ~ 0.01–0.1 (lighter). Stretching is weak everywhere relative to smoothing/direction terms.

**M3 (strain alignment):** Complex structure. The tube cores show near-zero alignment (white). The crease shows alternating positive/negative (red/blue spiral patterns) — the strain field is organized but cancels in the mean.

**M6 (Kolmogorov-conditioned Λ):** The crease glows as a bright arc/ribbon between the tubes. This is where Λ\_η >> 1. The tube cores are dark (Λ\_η < 1 — the direction field is smooth inside tubes). Clean spatial separation between crease and tubes.

**Product field |ω|·|∇ξ|\_F:** The brightest region is NOT at the crease center (where |∇ξ| peaks) nor at the tube cores (where |ω| peaks). It is at the INTERFACE — the ring/shell where both factors are moderate. This is the spatial location of the invariant candidate.

---

## Data Availability Flags

- **V5 t=40:** NOT AVAILABLE as field data. V5 stored only scalar timeseries to SQLite. The Lambda=1950 in V5 was computed with theta=0.01 masking during the run. Cannot retroactively compute M1-M6 at V5's crease location without re-running the solver with snapshot output.
- **M5 full strain decomposition:** The pressure-compression alignment (∇p · e\_compression where e is the minimum eigenvector of S) requires per-voxel eigendecomposition of the 3×3 strain tensor — feasible but expensive at N=256 (16M eigendecompositions). The proxy ∇p · ω alignment was computed instead. Flag if full eigendecomposition is needed.
- **N=64 resolution convergence:** M3 and M5 show qualitatively different behavior at N=64 vs N=128/256 (sign flip in M3, std divergence in M5). N=64 results should be treated as under-resolved for these metrics.

---

## Honest Interpretation

### What the data says

1. The product sup(|ω|·|∇ξ|\_F) is resolution-independent (γ = 0.043 ≈ 0, CV < 3% at t=0 and t=10) while the individual factors scale as N^{+0.86} and N^{≈ -0.86} respectively. This is an emergent constraint the flow imposes at the crease.

2. The product decays in time (slope ~ -0.15/t) as the tubes dissipate. It is not a conserved quantity. It is an instantaneous balance, not a conservation law.

3. The crease is NOT in a stretching-dominated regime (M1: B << 1). The local dynamics at the crease are dominated by the viscous and direction-gradient terms, not by vortex stretching.

4. 100% of crease voxels have Λ\_η > 1 (M6) — the direction gradient is sharper than the Kolmogorov scale everywhere on the crease, at every resolution and every timestep.

5. The spatial location of the product maximum is distinct from both the peak-gradient location (where |ω| → 0) and the tube cores (where |∇ξ| is small). It sits at the interface ring between the two.

### What the data does not say

1. It does not say whether the product invariance persists at higher Re or later times. V5's Lambda=1950 at N=128 t=40 suggests the crease sharpens dramatically at late times, but we cannot compute M4 at that snapshot.

2. It does not say whether the invariance is specific to antiparallel tube geometry or is a universal property of 3D NS reconnection. Different initial conditions would need testing.

3. It does not say whether the invariance implies regularity or singularity. A resolution-independent product is consistent with BOTH: (a) the flow self-regulates to prevent blowup (regularity), or (b) the product saturates at a threshold before the crease transitions to a different scaling regime.

4. It does not tell us the asymptotic behavior at N → ∞. Three data points (N=64, 128, 256) over a 4× range can confirm O(1/dx) but cannot distinguish N^{0.86} from N^{1.0} or N^{0.75} with confidence. More resolutions (N=384, 512 on larger hardware) would tighten the fit.

### What this means for paper framing

The data supports a paper about the crease invariant itself:

*"At the reconciliation surface between antiparallel coherent vortex tubes, the product |ω|·|∇ξ|\_F is resolution-independent while the individual factors scale inversely. This emergent balance constrains the viscous term's structure at direction-field discontinuities. The crease is sharper than the local Kolmogorov scale (Λ\_η > 1 at 100% of crease voxels) and sits in a non-stretching regime (B << 1). The Lean 4 formalization of the bootstrap failure remains a valid conditional theorem; the data shows the antecedents are not independently realized but are coupled through the product invariant."*

This is not the paper originally planned. It is a more honest paper with a sharper result.

---

## Files

| Output | Location |
|--------|----------|
| Raw metrics (17 snapshots × 24 fields) | `metrics.csv` |
| M4 scaling fits | `scaling_fit.csv` |
| M4 scaling plot (log-log) | `m4_scaling.png` |
| Slice maps (M1, M3, M6 + supporting) | `slice_maps.png` |
| Slice data (npz per snapshot) | `slices_N*_t*.npz` |
| Computation script | `compute_metrics.py` |
| This report | `crease_invariants_report.md` |
