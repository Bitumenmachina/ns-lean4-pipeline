# HDF5 Snapshots

The Phase 10b convergence snapshots (8.1GB total) are too large for the git repository.

## Files

| Snapshot | Size | Resolution | Time |
|----------|------|------------|------|
| run_64_t{0,10,20,30,36,40}.h5 | 23MB each | 64^3 | 6 snapshots |
| run_128_t{0,10,20,30,36,40}.h5 | 166-181MB each | 128^3 | 6 snapshots |
| run_256_t{0,10,20,25,28}.h5 | 1.3-1.5GB each | 256^3 | 5 snapshots |

Each HDF5 file contains 12 fields at the given resolution:
- `ux`, `uy`, `uz` — velocity components
- `omega_x`, `omega_y`, `omega_z` — vorticity components
- `omega_mag` — vorticity magnitude
- `xi_x`, `xi_y`, `xi_z` — direction field components
- `T1_norm`, `T2_norm` — decomposition term norms

XDMF metadata files are in `paraview/` — they reference HDF5 files by relative path.

## Obtaining the Data

Snapshots will be available as GitHub release artifacts or on Zenodo after publication.

To regenerate from scratch (~40 hours on 32GB consumer hardware):
```bash
cd solver/ffi-fftw/
maturin develop --release
cd python/
python3 phase10_convergence.py
```

## XDMF Files

The `.xdmf` metadata files live in `paraview/` (not this directory). Place
downloaded H5 files here and point `SNAPSHOT_DIR` at this directory when
running the ParaView scripts.
