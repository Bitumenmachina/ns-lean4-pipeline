# Reproducibility Guide

Step-by-step instructions to reproduce every result in this repository from a clean machine.

## Requirements

- Linux (tested on Fedora 43)
- 32 GB RAM minimum (for N=256 convergence runs)
- FFTW3 development libraries (`dnf install fftw3-devel` or `apt install libfftw3-dev`)
- Rust toolchain ≥ 1.94 (`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`)
- Python ≥ 3.12 with pip
- elan (Lean toolchain manager)

## Step 1: Install Dependencies (~5 min)

```bash
# Python
pip install numpy h5py matplotlib scipy z3-solver tabulate maturin

# Lean
curl -sSf -L https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y
source ~/.elan/env
```

## Step 2: Verify Lean Proof (~20 min first build, ~2 min cached)

```bash
cd lean/
lake build
# Expected: 964 jobs, 0 errors
grep -c '^\s*sorry' NSLean4/NSTheorems.lean
# Expected: 0
```

## Step 3: Verify Z3 Certificates (~1 sec)

```bash
cd ../z3/
python3 run_z3.py
# Expected: 10/10 SAT
```

## Step 4: Build Rust FFTW Solver (~2 min)

```bash
cd ../solver/ffi-fftw/
maturin develop --release
# Verify: python3 -c "import fftw_training; print('OK')"
```

## Step 5: Run Convergence Suite (~40 hours)

**This is the long step.** N=64 takes ~30 min, N=128 takes ~5 hours, N=256 takes ~35 hours.

```bash
cd python/
python3 phase10_convergence.py
# Output: run_64.csv, run_128.csv, run_256.csv + HDF5 snapshots in snapshots/
```

## Step 6: Run Crease Invariant Analysis (~30 min)

Requires HDF5 snapshots from Step 5.

```bash
cd ../../../data/crease_invariants/
SNAPSHOT_DIR=/path/to/snapshots python3 compute_metrics.py
# Output: metrics.csv, scaling_fit.csv, slice maps
```

## Step 7: Verify M4 Result

```bash
python3 -c "
import csv
rows = list(csv.DictReader(open('scaling_fit.csv')))
for r in rows:
    print(f't={r[\"t\"]:>3}  β={float(r[\"beta\"]):>6.4f}  γ={float(r[\"product_slope\"]):>7.4f}  R²={float(r[\"R2_beta\"]):>6.4f}')
beta_mean = sum(float(r['beta']) for r in rows) / len(rows)
gamma_mean = sum(float(r['product_slope']) for r in rows) / len(rows)
print(f'β mean = {beta_mean:.4f}  (expect ~0.86)')
print(f'γ mean = {gamma_mean:.4f}  (expect ~0.04, ≈ 0)')
"
```

## Step 8: Generate Visualizations

```bash
# Convergence overlay (requires matplotlib)
cd ../../  # back to repo root
# ... convergence plot generation is embedded in the convergence runner

# ParaView crease visualization (requires ParaView + HDF5 snapshots)
cd paraview/
SNAPSHOT_DIR=/path/to/snapshots pvbatch convergence_multires.py --time 20
# Output: output/convergence_t20.png
```

## Time Budget (Consumer Hardware: i7-14700HX, 32GB, RTX 4060)

| Step | Time | Notes |
|------|------|-------|
| Dependencies | 5 min | One-time |
| Lean build | 20 min | First build downloads Mathlib; cached ~2 min after |
| Z3 certificates | <1 sec | |
| Rust solver build | 2 min | |
| N=64 convergence | 30 min | |
| N=128 convergence | 5 hours | |
| N=256 convergence | 35 hours | May thermal-throttle on laptops |
| Crease invariant analysis | 30 min | Requires all snapshots |
| **Total** | **~42 hours** | Dominated by N=256 run |

## Pinned Versions

| Component | Version | Pin |
|-----------|---------|-----|
| Lean 4 | v4.29.0 | `lean-toolchain` |
| Mathlib | v4.29.0 | `lake-manifest.json` (commit `8a178386...`) |
| Z3 | 4.16.0 | `pip install z3-solver==4.16.0` |
| Rust | 1.94.1 | `rust-toolchain.toml` |
| FFTW | 3.3.10 | System package |
| Python | ≥3.12 | `pyproject.toml` |
| NumPy | ≥2.0 | `requirements.txt` |

## Not in This Repo

`pipeline/extract_params.py` reads from source SQLite databases (`~/ns_archive/databases/`) that are not distributed. The extracted parameters are already committed as `params.json` and the Z3/Lean pipeline runs from that file. You only need the source databases to re-extract from raw simulation output.
