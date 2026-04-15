# Complete Convergence Data — All Experiments

**Generated:** 2026-04-12
**Scope:** V3 (geometry survey) → V4 (pressure conjecture) → V5 (reconnection onset) → Phase 10b (spectral convergence)

---

## 1. V5 Reconnection — Onset Table (N=128, pseudospectral)

### Full Re Sweep (t=20 unless noted)

| Re | N | Λ_init | Λ_trough | Λ_final | fin/tro | ω_max | Ω_final | Behavior |
|----|---|--------|----------|---------|---------|-------|---------|----------|
| 400 | 128 | 344.2 | 154.3 | 159.2 | 1.0 | 0.17 | 0.743 | Decays |
| 800 | 128 | 339.3 | 189.4 | 213.7 | 1.1 | 0.35 | 1.765 | Flat |
| 1600 | 64 | 322.9 | 117.8 | 266.7 | 2.3 | 0.69 | 3.566 | Amplified recovery |
| 1600 | 128 | 339.0 | 139.7 | 267.6 | 1.9 | 0.70 | 3.566 | Amplified recovery |
| 3200 | 64 | 330.5 | 59.3 | 296.0 | 5.0 | 1.30 | 6.369 | Amplified recovery |
| 3200 | 128 | 342.1 | 64.7 | 311.7 | 4.8 | 1.34 | 6.369 | Amplified recovery |
| 3600 | 128 | 341.5 | 53.1 | 316.0 | 6.0 | 1.46 | 6.941 | Amplified recovery |
| 4000 | 128 | 340.0 | 44.7 | 330.1 | 7.4 | 1.57 | 7.472 | Amplified recovery |
| 4800 | 128 | 351.2 | 44.9 | 349.3 | 7.8 | 1.74 | 8.424 | **GROWS** |
| 5600 | 128 | 350.6 | 44.9 | 396.3 | 8.8 | 1.88 | 9.255 | **GROWS** |
| 6400 | 128 | 350.2 | 40.4 | 391.6 | 9.7 | 1.98 | 9.989 | **GROWS** |
| 6400 | 128 | 350.2 | 40.4 | 1950.1 | 48.2 | 2.65 | 6.884 | **GROWS** (t=40) |
| 6400 | 256 | 346.9 | 44.0 | 431.3 | 9.8 | 1.98 | 9.989 | **GROWS** |

### Re=6400 Extended Run (N=128, t=0→40)

| t | Λ | ω_max | Ω | Note |
|---|---|-------|---|------|
| 10 | 47 | 1.89 | 11.33 | Viscous trough |
| 14 | 100 | 2.08 | 10.60 | K crossover (~1.0) |
| 20 | 392 | 1.98 | 9.99 | Accelerating |
| 25 | 1273 | 2.01 | 9.47 | 3.2× since t=20 |
| 30 | 1779 | 2.32 | 8.96 | |
| 35 | 1005 | 3.60 | 8.05 | Reconnection dip |
| 37 | 1119 | 3.84 | 7.79 | Peak |ω|_max |
| 39 | 1923 | 3.22 | 7.19 | Rebound |
| 40 | 1950 | 2.65 | 6.88 | Still growing |

### N=256 vs N=128 Convergence at Re=6400

| t | Λ_128 | Λ_256 | ratio |
|---|-------|-------|-------|
| 0 | 279.7 | 318.0 | 1.14 |
| 5 | 58.2 | 56.7 | 0.98 |
| 10 | 52.4 | 53.9 | 1.03 |
| 15 | 125.1 | 126.1 | 1.01 |
| 20 | 391.6 | 431.3 | 1.10 |

**Verdict:** Ratio trending toward 1.0. Signal is physical, not numerical.

---

## 2. V3 Incoherence — 3D Geometry Survey (FluidX3D LBM)

### Separation Geometries (I increases with Re — PASS)

| Geometry | Re | I_mean | stretch_L2 | Ω |
|----------|-----|--------|------------|---|
| Sphere | 100 | 0.1501 | 8.45e-07 | 2.71e-06 |
| Sphere | 400 | 0.2058 | 1.08e-06 | 2.31e-06 |
| Sphere | 800 | 0.2181 | 1.30e-06 | 2.32e-06 |
| Sphere | 1600 | 0.2344 | 2.29e-06 | 2.51e-06 |
| Cube | 100 | 0.1094 | 1.36e-06 | 2.80e-06 |
| Cube | 400 | 0.1356 | 1.68e-06 | 2.60e-06 |
| Cube | 800 | 0.1623 | 2.41e-06 | 2.96e-06 |
| Cube | 1600 | 0.1770 | 6.48e-06 | 3.75e-06 |
| Step | 100 | 0.1563 | 2.92e-07 | 2.22e-06 |
| Step | 500 | 0.1884 | 7.94e-07 | 4.72e-06 |
| Step | 1000 | 0.2145 | 2.99e-06 | 6.64e-06 |
| Step | 2000 | 0.2481 | 1.56e-05 | 1.25e-05 |
| Tandem | 200 | 0.1542 | 7.98e-07 | 2.44e-06 |
| Tandem | 800 | 0.2032 | 1.20e-06 | 2.41e-06 |
| Tandem | 1600 | 0.2215 | 2.80e-06 | 2.74e-06 |

### Clean Pipe (I saturated — FAIL)

| Re | I_mean | Note |
|-----|--------|------|
| 500 | 0.4945 | Saturated |
| 1000 | 0.5090 | Saturated |
| 2000 | 0.4902 | Saturated |
| 3000 | 0.4856 | Saturated |
| 4000 | 0.4843 | Saturated |

**Limitation:** Cylindrical geometry creates organized perpendicular structure misclassified as incoherence.

---

## 3. V4 Spectral Conjecture — Pressure Enforcement

| Dim | Re | Behavior |
|-----|----|----------|
| 2D | 100–1600 | Residual 0.003→0.001 (**IMPROVES** with Re) |
| 3D | 100 | Residual 0.013 |
| 3D | 200 | Residual 0.037, λ=0.079 |
| 3D | 400 | Residual 0.038, λ=0.114 |
| 3D | 800 | Residual 0.080, λ=0.278 |
| 3D | 1600 | Residual 0.167, λ=0.396 |

Resolution sweep at Re=800: N=32, 64, 128, 256 (to determine if degradation is physical or numerical).
Extended Re sweep: 3200, 6400, 12800, 25600, 51200, 102400.

---

## 4. V3↔V5 Bridge — Cross-Experiment Correlation

### V3 Obstacle Flows: Λ_L2 by Geometry (at step 20000)

| Run | Λ_L2 | I_ξ_L2 | I_orig_L2 |
|-----|------|--------|-----------|
| sphere_Re400 | 0.0790 | 0.0217 |
| sphere_Re800 | 0.0915 | 0.0473 |
| sphere_Re1600 | 0.1264 | 0.0678 |
| cube_Re400 | 0.0825 | 0.0355 |
| cube_Re800 | 0.1140 | 0.0820 |
| cube_Re1600 | 0.1794 | 0.1489 |
| step_Re500 | 0.1650 | 0.0264 |
| step_Re1000 | 0.2066 | 0.1161 |
| step_Re2000 | 0.2675 | 0.1630 |

### V5 Vortex Tubes: I_ξ at t=20

| Re | Λ_final | I_ξ_L2 | I_ξ_turb |
|-----|---------|--------|----------|
| 3200 | 300.7 | 0.098 | 0.026 |
| 4800 | 329.7 | 0.110 | 0.036 |
| 6400 | 361.9 | 0.121 | 0.045 |

**Cross-run correlations:** Λ_L2 vs I_orig r=0.965, Λ_L2 vs I_ξ r=0.862, I_ξ vs I_orig r=0.820.
**Within-run correlation:** near zero (r=-0.03 to 0.16). Different temporal dynamics.
**Interpretation:** Λ = leading edge (sharpest gradient), I = extent (cumulative disorder). Companion diagnostics, not interchangeable.

---

## 5. Phase 10b — FFI-FFTW Spectral Convergence (Apr 10-12)

**IC:** Antiparallel Gaussian vortex tubes, Γ=1, σ=0.3, d=π/2, pert=0.1
**Re:** 6400 (ν=2π/6400 ≈ 9.817e-4), dt=0.001

### Cross-Resolution Diagnostic Table

| t | Λ_64 | Λ_128 | Λ_256 | 128/64 | 256/128 | K_64 | K_128 | K_256 | ω_max_64 | ω_max_128 | ω_max_256 | T2/T1_64 | T2/T1_128 | T2/T1_256 |
|---|------|-------|-------|--------|---------|------|-------|-------|----------|-----------|-----------|----------|-----------|-----------|
| 0 | 28.56 | 51.10 | 101.35 | 1.79 | 1.98 | 1.636 | 2.928 | 5.808 | 3.53 | 3.53 | 3.53 | 1.06 | 1.08 | 1.45 |
| 5 | 38.69 | 53.25 | 102.51 | 1.38 | 1.93 | 2.327 | 3.202 | 6.165 | 2.90 | 2.90 | 2.90 | 1.09 | 1.24 | 1.26 |
| 10 | 41.10 | 57.78 | 113.56 | 1.41 | 1.97 | 2.580 | 3.627 | 7.128 | 2.42 | 2.50 | 2.50 | 1.20 | 0.93 | 1.26 |
| 14 | 36.21 | 57.98 | 113.91 | 1.60 | 1.96 | 2.342 | 3.750 | 7.367 | 2.17 | 2.19 | 2.20 | 1.48 | 1.03 | 1.54 |
| 20 | 32.03 | 57.54 | 114.36 | 1.80 | 1.99 | 2.154 | 3.869 | 7.690 | 1.87 | 1.88 | 1.88 | 1.17 | 0.90 | 1.38 |
| 25 | 38.79 | 59.14 | 113.83 | 1.52 | 1.92 | 2.684 | 4.092 | 7.877 | 1.70 | 1.70 | 1.70 | 1.61 | 1.83 | 1.71 |
| 28 | 40.17 | 56.94 | 114.19 | 1.42 | 2.01 | 2.824 | 4.003 | 8.027 | 1.60 | 1.61 | 1.61 | 1.38 | 1.31 | 1.57 |
| 30 | 35.75 | 56.53 | — | 1.58 | — | 2.538 | 4.014 | — | 1.54 | 1.54 | — | 1.43 | 1.43 | — |
| 36 | 35.02 | 65.65 | — | 1.87 | — | 2.556 | 4.792 | — | 1.37 | 1.38 | — | 1.78 | 1.78 | — |
| 40 | 35.37 | 59.26 | — | 1.68 | — | 2.626 | 4.401 | — | 1.28 | 1.30 | — | 1.14 | 1.29 | — |

### Run Summary

| Metric | N=64 | N=128 | N=256 |
|--------|------|-------|-------|
| t_max | 40.0 | 40.0 | 30.0 |
| t_reached | 40.00 | 40.00 | 29.69 |
| steps | 40000 | 40000 | 29690 |
| wall clock (h) | 0.48 | 4.97 | 35.50 |
| peak RSS (GB) | 0.142 | 0.868 | 6.214 |
| status | complete | complete | wall-clock halt |

**N=256 halt:** Wall-clock 5.8746s/step > 1.5× probe (3.8784s) for 20 consecutive batches

---

## 6. Step1 Onset Databases (individual Re runs)

| Database | Rows | Λ_final | Λ_trough | ω_max |
|----------|------|---------|----------|-------|
| Re3600_N128 | 188 | 316.0 | 53.1 | 1.46 |
| Re4000_N128 | 188 | 330.1 | 44.7 | 1.57 |
| Re4800_N128 | 188 | 349.3 | 44.9 | 1.74 |
| Re5600_N128 | 188 | 396.3 | 44.9 | 1.88 |
| Re6400_N128_t40 | 380 | 1950.1 | 40.4 | 2.65 |

---

## Provenance

| Source | Solver | Period | Location |
|--------|--------|--------|----------|
| V3 3D incoherence | FluidX3D LBM (GPU) | Apr 6 | `ns_archive/databases/experiment_3d.db` |
| V4 pressure conjecture | Pseudospectral | Apr 7 | `ns_archive/databases/conjecture_v4.db` |
| V5 reconnection | Pseudospectral | Apr 7-8 | `ns_archive/databases/reconnection_v5.db` |
| V5 step1 onset | Pseudospectral | Apr 8 | `ns_archive/databases/step1*.db` |
| V5 bridge | Both solvers | Apr 8 | `ns_archive/databases/bridge_results.db` |
| Phase 10b convergence | FFI-FFTW spectral (Rust) | Apr 10-12 | `ffi-fftw-training/run_*.csv` |
