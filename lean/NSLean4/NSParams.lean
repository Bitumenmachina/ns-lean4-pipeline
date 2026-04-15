import Mathlib.Data.Rat.Defs
import Mathlib.Algebra.Order.Field.Rat

/-!
  NSParams.lean — Auto-generated from V5 experimental data
  Source: params.json (6487 bytes)
  Databases: experiment.db, experiment_v2.db, bridge_results.db, experiment_3d.db, conjecture_v4.db, reconnection_v5.db, step1_3D_Re4800_N128.db, step1_3D_Re5600_N128.db, step1_3D_Re6400_N128_t40.db, step1b_3D_Re3600_N128.db, step1b_3D_Re4000_N128.db

  K = Λ · η  where  η = (ν³/ε)^(1/4)  (Kolmogorov scale)
  All values are exact rationals via limit_denominator(10^8).
-/

namespace NS

-- Solver parameters
noncomputable def nu_Re6400 : ℚ := (3927 / 4000000 : ℚ)

-- Key measurements at Re=6400
def Lambda_Re6400_t20 : ℚ := (39159 / 100 : ℚ)
def Lambda_Re6400_t40 : ℚ := (195007 / 100 : ℚ)
def Lambda_Re6400_N256_t20 : ℚ := (21567 / 50 : ℚ)
def eta_Re6400_t40 : ℚ := (1209 / 62500 : ℚ)
def epsilon_Re6400_t40 : ℚ := (3379 / 500000 : ℚ)
def K_Re6400_t40 : ℚ := (943 / 25 : ℚ)
def K_Re6400_t20 : ℚ := (69 / 10 : ℚ)

-- Enstrophy (monotonically decreasing)
def Omega_Re6400_t0 : ℚ := (1231 / 100 : ℚ)
def Omega_Re6400_t20 : ℚ := (999 / 100 : ℚ)
def Omega_Re6400_t40 : ℚ := (68841 / 10000 : ℚ)

-- Resolution convergence ratios (must be < 1.3 and trending → 1)
def conv_ratio_128_64 : ℚ := (11701 / 10000 : ℚ)
def conv_ratio_256_128 : ℚ := (2199 / 2000 : ℚ)

-- Onset table: (Re, Lambda_final, final/trough)
-- Monotonically increasing final/trough with Re
def onset_table : List (ℚ × ℚ × ℚ) :=
  [((400 : ℚ), (159 : ℚ), (9 / 10 : ℚ)),
  ((800 : ℚ), (214 : ℚ), (11 / 10 : ℚ)),
  ((1600 : ℚ), (268 : ℚ), (19 / 10 : ℚ)),
  ((3200 : ℚ), (312 : ℚ), (24 / 5 : ℚ)),
  ((3600 : ℚ), (316 : ℚ), (6 : ℚ)),
  ((4000 : ℚ), (330 : ℚ), (7 : ℚ)),
  ((4800 : ℚ), (349 : ℚ), (39 / 5 : ℚ)),
  ((5600 : ℚ), (396 : ℚ), (44 / 5 : ℚ)),
  ((6400 : ℚ), (392 : ℚ), (49 / 5 : ℚ))]

-- Staircase growth at Re=6400: (t, Lambda, K=Lambda*eta)
-- Pattern: oscillatory ratchet, each peak higher
def staircase : List (ℚ × ℚ × ℚ) :=
  [((239 / 50 : ℚ), (5421 / 100 : ℚ), (23 / 25 : ℚ)),
  ((261 / 50 : ℚ), (1278 / 25 : ℚ), (87 / 100 : ℚ)),
  ((268 / 25 : ℚ), (4741 / 100 : ℚ), (81 / 100 : ℚ)),
  ((701 / 50 : ℚ), (10039 / 100 : ℚ), (87 / 50 : ℚ)),
  ((20 : ℚ), (39159 / 100 : ℚ), (69 / 10 : ℚ)),
  ((126 / 5 : ℚ), (127343 / 100 : ℚ), (91 / 4 : ℚ)),
  ((30 : ℚ), (177887 / 100 : ℚ), (3221 / 100 : ℚ)),
  ((888 / 25 : ℚ), (25126 / 25 : ℚ), (187 / 10 : ℚ)),
  ((1831 / 50 : ℚ), (22381 / 20 : ℚ), (2099 / 100 : ℚ)),
  ((1941 / 50 : ℚ), (9614 / 5 : ℚ), (3679 / 100 : ℚ)),
  ((998 / 25 : ℚ), (195007 / 100 : ℚ), (943 / 25 : ℚ))]

-- K = Λ·η at t=20 for each Re (K decreases with Re at fixed t)
def K_vs_Re : List (ℚ × ℚ) :=
  [((400 : ℚ), (43 / 2 : ℚ)),
  ((800 : ℚ), (82 / 5 : ℚ)),
  ((1600 : ℚ), (61 / 5 : ℚ)),
  ((3200 : ℚ), (87 / 10 : ℚ)),
  ((3600 : ℚ), (81 / 10 : ℚ)),
  ((4000 : ℚ), (79 / 10 : ℚ)),
  ((4800 : ℚ), (37 / 5 : ℚ)),
  ((5600 : ℚ), (38 / 5 : ℚ)),
  ((6400 : ℚ), (69 / 10 : ℚ))]

-- Bifurcation regime boundaries (continuous, not sharp)
def Re_regime1_upper : ℚ := (800 : ℚ)
def Re_regime2_lower : ℚ := (1600 : ℚ)
def Re_regime2_upper : ℚ := (4800 : ℚ)
def Re_regime3_lower : ℚ := (5600 : ℚ)

/-- The kill metric K = Λ · η. Dimensionless.
    Direction field gradient (Λ) scaled by Kolmogorov length (η).
    K > 1 means direction varies faster than dissipation scale. -/
noncomputable def K_metric (Lambda eta : ℚ) : ℚ := Lambda * eta

end NS
