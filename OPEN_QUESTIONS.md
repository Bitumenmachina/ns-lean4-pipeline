# Open Questions

## 1. Is the crease invariant a barrier to singularity or a signature of one?

The product |ω|·sup|∇ξ|\_F is resolution-independent (γ = 0.04 ≈ 0) at the interface ring between antiparallel vortex tubes at Re=6400. Two interpretations are consistent with this data:

**Interpretation A (regularity barrier):** The flow self-regulates — wherever |∇ξ| grows, |ω| decays to compensate, preventing the viscous term from becoming distributional. The invariant is a manifestation of the smoothing mechanism succeeding at the crease.

**Interpretation B (singularity precursor):** The invariant is the balanced state that holds during the approach to reconnection. At sufficiently late time or high Re, the balance breaks — one factor grows faster than the other can compensate — and the product diverges. The invariant is a plateau before a transition.

**The authors' current hypothesis is Interpretation B**, based on the V5 observation that Λ accelerates through reconnection-relaxation cycles at Re=6400 t=40 (Λ = 1950 at N=128, still growing). However, V5 t=40 field data is not available for the M4 analysis — V5 stored only scalar timeseries. Resolving this requires:

- Re-running V5 with field snapshot output at t=30, 35, 40 (N=128)
- Running Phase 10b to t=40 at N=256 (requires ~48h continuous compute; current hardware thermally throttles at ~35h)
- Ideally: N=256 on hardware with active cooling or N=512 on 64GB+ RAM

**Status:** Deferred pending access to larger compute. The paper should present the invariant as observed and explicitly state this question as open.

## 2. Is the invariant specific to antiparallel tube reconnection?

All current data uses antiparallel Gaussian vortex tubes (Γ=1, σ=0.3, d=π/2). The invariant could be:
- **Universal:** a property of all 3D NS reconnection events
- **Geometry-dependent:** specific to antiparallel tubes by symmetry
- **IC-dependent:** an artifact of the smooth Gaussian initial condition

Testing requires different ICs: orthogonal tubes, trefoil knots, random vorticity fields, Taylor-Green vortex at higher Re.

## 3. Does the invariant persist at higher Re?

Current data is at Re=6400 only. The V5 onset table shows the crease intensifies continuously from Re=1600 through Re=6400 with no sharp transition. Whether the product invariant holds at Re=12800 or Re=25600 (where V4 ran but without direction-field diagnostics) is unknown.

## 4. What is the precise exponent β?

β = 0.860 ± 0.051 from three resolutions (N=64, 128, 256). Three data points cannot distinguish β=0.86 from β=1.0 (exact O(1/dx)) or β=0.75. N=384 and N=512 would tighten the fit. β=1 has a clean interpretation (pointwise discontinuity in ξ); β<1 suggests the crease has finite width that scales sub-linearly.
