# ns-lean4-pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19584003.svg)](https://doi.org/10.5281/zenodo.19584003)

An end-to-end verification pipeline: DNS experiment data &rarr; exact rational parameters &rarr; Z3 certificates &rarr; Lean 4 proofs. Zero `sorry`, zero manual transcription, zero trust gaps between simulation and formal mathematics.

The pipeline is solver-agnostic. This repository includes a conditional viscosity bootstrap failure theorem as a worked example.

**What this is not.** Not a Clay Millennium resolution. Not a blowup proof. Not unconditional anything. A reproducible geometric observation, formally verified where it can be, honestly limited where it cannot.

## Pipeline Architecture

```
SQLite DBs ──► extract_params.py ──► params.json
                                         │
                          ┌───────────────┼───────────────┐
                          ▼               ▼               ▼
                   s1_gen_lean_defs.py  s7_z3_certs.py  (human review)
                          │               │
                          ▼               ▼
                   NSParams.lean    all_certs.smt2
                   (21 exact          (10 QF_LRA
                    rationals)         blocks)
                                        │
                                        ▼
                                   run_z3.py
                                        │
                                        ▼
                                  results.json
                                   (10/10 SAT)
                                        │
                                        ▼
                                  s8_z3_to_lean.py
                                        │
                                        ▼
                                  NSCerts.lean
                                  (10 axioms)
                                        │
                          ┌─────────────┘
                          ▼
                   NSTheorems.lean
                   (13 theorems,
                    0 sorry)
                          │
                          ▼
                     lake build
                   (964 jobs, 0 errors)
```

Every number in the formal proof traces back to a specific SQL query on experimental data. Every axiom traces to a Z3-certified arithmetic bound. The trust model is documented in [AUDIT.md](AUDIT.md).

## Worked Example: Navier-Stokes Direction-Field Crease Invariant

Decomposing vorticity as &omega; = |&omega;|&xi;, the viscous term &nu;&nabla;&sup2;&omega; contains |&omega;|&nabla;&sup2;&xi;. At the reconciliation surface between antiparallel vortex tubes, the direction field &xi; = &omega;/|&omega;| develops a crease where Lip(&xi;) scales as O(1/dx).

**Conditional theorem (Lean 4, 0 sorry):** IF &Lambda; = Lip(&xi;) diverges AND |&omega;| > 0 at the crease, THEN &nu;&nabla;&sup2;&omega; is unbounded in L&sup2;.

**Empirical finding (M4 scaling invariant):** At the crease interface, the product sup(|&omega;|&middot;|&nabla;&xi;|) is resolution-independent (&gamma; = 0.04 &asymp; 0) while the individual factors scale as N^{+0.86} and N^{-0.86} respectively. The antecedents H1 and H2 are coupled through this product invariant and are not independently realized in the current simulations. See [data/crease_invariants/](data/crease_invariants/) for the full analysis.

Key observations from pseudospectral DNS at Re=6400:
- &Lambda; = sup|&nabla;&xi;|\_F scales as N^{0.86 &plusmn; 0.05} (R&sup2; > 0.96)
- Product sup(|&omega;|&middot;|&nabla;&xi;|) ~ N^{0.04} (CV < 3% at t=0,10)
- 100% of crease voxels have &Lambda;&eta; > 1 (direction gradient sharper than Kolmogorov scale)
- Crease is in non-stretching regime (balance ratio B &lt;&lt; 1)
- See [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) for unresolved interpretive questions

## Repository Structure

```
lean/               Lean 4 formalization (lakefile.toml, Mathlib v4.29.0)
z3/                 Z3 certificates (10 QF_LRA blocks, 10/10 SAT)
solver/             DNS solvers (V5 pseudospectral + Rust FFTW spectral)
pipeline/           Data extraction and code generation scripts
data/
  convergence/      Phase 10b cross-resolution tables and overlay plot
  crease_invariants/  M1–M6 metric outputs, scaling fits, slice maps
  onset/            V5 9-Re sweep summary
  bridge/           V3↔V5 cross-experiment correlation data
  snapshots/        README pointing to Zenodo/release for H5 (>50MB)
paraview/           ParaView visualization scripts for crease rendering
paper/              LaTeX source + CLAUDE_NOTES.md
AUDIT.md            Axiom-to-certificate trail with trust boundaries
NOTATION.md         Symbol definitions across paper, code, and data
OPEN_QUESTIONS.md   Unresolved interpretive questions (barrier vs. signature)
REPRODUCE.md        Step-by-step reproducibility guide with time estimates
```

## Reproduce Everything

**Requirements:** Linux, 32GB RAM, FFTW3, Rust toolchain, Python &ge;3.12

```bash
# 1. Verify the formal proof (~20 min on first build, downloads Mathlib)
cd lean/
elan install                     # Lean toolchain
lake build                       # 964 jobs, 0 errors, 0 sorry

# 2. Verify Z3 certificates (~1 sec)
pip install z3-solver
python3 ../z3/run_z3.py          # 10/10 SAT

# 3. Regenerate everything from source databases (optional)
cd ../pipeline/
python3 extract_params.py        # SQLite → params.json
python3 s1_gen_lean_defs.py      # params.json → NSParams.lean
python3 s7_z3_certs.py           # params.json → all_certs.smt2
python3 ../z3/run_z3.py          # certificates → results.json
python3 s8_z3_to_lean.py         # results.json → NSCerts.lean
cd ../lean/ && lake build        # rebuild with regenerated files

# 4. Run the DNS solver (~40 hours for full convergence suite)
cd ../solver/ffi-fftw/
maturin develop --release        # build Rust FFTW kernel
cd python/
python3 phase10_convergence.py   # N=64/128/256 convergence runs

# 5. Visualize convergence (requires ParaView)
cd ../../paraview/
pvbatch convergence_multires.py  # generates output/convergence_t20.png
```

**Consumer hardware tested on:** Lenovo Legion 7 (i7-14700HX, 32GB RAM, RTX 4060). Full convergence suite completes in ~40 hours. Lean build completes in ~20 minutes. Z3 verification completes in <1 second.

## Pinned Versions

| Component | Version | Commit |
|-----------|---------|--------|
| Lean 4 | v4.29.0 | `leanprover/lean4:v4.29.0` |
| Mathlib | v4.29.0 | `8a178386ffc0f5fef0b77738bb5449d50efeea95` |
| Z3 (Python) | 4.16.0 | |
| Rust | 1.94.1 | |
| FFTW | 3.3.10 | |
| Python | &ge;3.12 | |

## Trust Model

| Layer | Verified By | Status |
|-------|------------|--------|
| Algebraic implications (13 theorems) | Lean 4 kernel | 0 sorry, 0 native_decide |
| Arithmetic bounds (10 certificates) | Z3 QF_LRA (decidable, sound) | 10/10 SAT |
| Empirical antecedents (H1: &Lambda;&rarr;&infin;, H2: \|&omega;\|>0) | Pseudospectral DNS | Supported, not provable from finite data |

See [AUDIT.md](AUDIT.md) for the complete axiom-to-certificate mapping.

## Citation

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

[MIT](LICENSE)
