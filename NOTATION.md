# Notation

Consistent notation across paper, code, and data files.

## Fields

| Symbol | Code name | Description |
|--------|-----------|-------------|
| **u** | `ux, uy, uz` | Velocity field |
| **ω** = ∇ × u | `omega_x, omega_y, omega_z` | Vorticity |
| \|ω\| | `omega_mag` | Vorticity magnitude |
| **ξ** = ω/\|ω\| | `xi_x, xi_y, xi_z` | Vorticity direction field |
| S = ½(∇u + ∇uᵀ) | computed | Rate-of-strain tensor |
| p | computed from Poisson | Pressure (zero-mean gauge) |

## Crease Diagnostics

| Symbol | Code name | Definition |
|--------|-----------|------------|
| Λ | `Lambda`, `grad_xi_F` | sup\|∇ξ\|\_F (Lipschitz constant of direction field) |
| K | `K` | Λ · η (crease sharpness relative to Kolmogorov scale) |
| η | `eta_local` | (ν³/ε)^{1/4} (Kolmogorov microscale) |
| ε | `epsilon_local` | 2ν S:S (local dissipation rate) |
| Λ\_η | `Lambda_eta` | \|∇ξ\| · η (pointwise Kolmogorov-conditioned Λ) |
| T1 | `T1_norm` | \|∇\|ω\|\| (magnitude gradient term) |
| T2 | `T2_norm` | \|ω\| · \|∇²ξ\| (direction gradient term) |

## Scaling Exponents

| Symbol | Meaning | Measured Value |
|--------|---------|----------------|
| β | sup\|∇ξ\|\_F ~ N^β | 0.860 ± 0.051 |
| α (implied) | \|ω\| at product-max location ~ N^α | ≈ −β |
| γ = α + β | sup(\|ω\|·\|∇ξ\|) ~ N^γ | 0.043 ≈ 0 |

## Metric Labels (M1–M6)

| Label | Quantity | Purpose |
|-------|----------|---------|
| M1 | B = \|ω·∇u\| / (ν\|∇²ω\| + \|ω\|·\|∇ξ\|) | Local balance ratio |
| M2 | h = u · ω | Helicity density |
| M3 | cos(θ) = (S·ω)·ω / (\|S·ω\|·\|ω\|) | Strain-vorticity alignment |
| M4 | \|ω\| · \|∇ξ\|\_F at product maximum | **Invariant candidate** |
| M5 | ∇p and its alignment | Pressure gradient at crease |
| M6 | Λ\_η = \|∇ξ\| · η | Kolmogorov-conditioned direction gradient |

## Physical Parameters

| Symbol | Value | Definition |
|--------|-------|------------|
| Re | 6400 | Reynolds number |
| ν | 2π/6400 ≈ 9.817 × 10⁻⁴ | Kinematic viscosity |
| L | 2π | Domain size (periodic cube) |
| dt | 0.001 | Timestep |
| Γ | 1 | Initial tube circulation |
| σ | 0.3 | Initial tube core radius |
| d | π/2 | Initial tube separation |
