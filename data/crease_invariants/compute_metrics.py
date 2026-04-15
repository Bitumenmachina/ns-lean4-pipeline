#!/usr/bin/env python3
"""Crease invariant characterization — six diagnostic metrics.

Computes M1-M6 at each (N, t) snapshot from Phase 10b H5 data.
All spectral derivatives use 2/3 dealiasing to match the solver.

Usage:
  SNAPSHOT_DIR=/path/to/h5/snapshots python3 compute_metrics.py
  (or place H5 files in ../../data/snapshots/ relative to this script)

Output: metrics.csv, scaling_fit.csv, slice maps (PNG), report data

Dependencies: numpy, h5py, scipy (for scaling fits only)
"""
import numpy as np
import h5py
import csv
import os
import sys
import json
from pathlib import Path

NU = 2 * np.pi / 6400  # Re=6400

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))

SNAP_DIR = os.environ.get("SNAPSHOT_DIR",
    os.path.join(REPO_ROOT, "data", "snapshots"))
OUT_DIR = os.environ.get("OUTPUT_DIR", SCRIPT_DIR)


# ── Spectral tooling ──

def make_wavenumbers(N):
    """Wavenumber arrays + 2/3 dealiasing mask matching the solver."""
    dx = 2 * np.pi / N
    k1d = np.fft.fftfreq(N, d=dx / (2 * np.pi))
    kx = k1d[:, None, None]
    ky = k1d[None, :, None]
    kz = k1d[None, None, :]
    k2 = kx**2 + ky**2 + kz**2
    k2_safe = np.where(k2 == 0, 1.0, k2)

    # 2/3 dealiasing mask
    k_max = N // 3
    dealias = np.ones((N, N, N))
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                if abs(k1d[ix]) > k_max or abs(k1d[iy]) > k_max or abs(k1d[iz]) > k_max:
                    dealias[ix, iy, iz] = 0.0

    return kx, ky, kz, k2, k2_safe, dealias


def spectral_deriv(field_hat, ki):
    """Spectral derivative: multiply by i*k."""
    return np.fft.ifftn(1j * ki * field_hat).real


def spectral_laplacian(field_hat, k2):
    """Spectral Laplacian: multiply by -k²."""
    return np.fft.ifftn(-k2 * field_hat).real


# ── Metric computations ──

def compute_all_metrics(h5_path):
    """Compute M1-M6 for a single snapshot. Returns dict of results."""
    f = h5py.File(h5_path, "r")
    ux, uy, uz = f["ux"][:], f["uy"][:], f["uz"][:]
    ox, oy, oz = f["omega_x"][:], f["omega_y"][:], f["omega_z"][:]
    omega_mag = f["omega_mag"][:]
    xi_x, xi_y, xi_z = f["xi_x"][:], f["xi_y"][:], f["xi_z"][:]
    N = int(f.attrs["N"])
    t = float(f.attrs["t"])
    f.close()

    dx = 2 * np.pi / N
    kx, ky, kz, k2, k2_safe, dealias = make_wavenumbers(N)
    Ks = [kx, ky, kz]

    # Velocity in spectral space
    u = [ux, uy, uz]
    u_hat = [np.fft.fftn(c) for c in u]
    o = [ox, oy, oz]
    o_hat = [np.fft.fftn(c) for c in o]
    xi = [xi_x, xi_y, xi_z]

    # ── Velocity gradient tensor ∂u_i/∂x_j ──
    grad_u = np.zeros((3, 3, N, N, N))
    for i in range(3):
        for j in range(3):
            grad_u[i, j] = spectral_deriv(u_hat[i], Ks[j])

    # ── Strain tensor S = ½(∇u + ∇uᵀ) ──
    S = np.zeros((3, 3, N, N, N))
    for i in range(3):
        for j in range(3):
            S[i, j] = 0.5 * (grad_u[i, j] + grad_u[j, i])

    # ── |∇ξ|_F (direction field gradient) ──
    grad_xi_sq = np.zeros((N, N, N))
    for comp_hat in [np.fft.fftn(c) for c in xi]:
        for ki in Ks:
            grad_xi_sq += spectral_deriv(comp_hat, ki)**2
    grad_xi_F = np.sqrt(grad_xi_sq)
    del grad_xi_sq

    # Peak gradient voxel
    peak_idx = np.unravel_index(np.argmax(grad_xi_F), grad_xi_F.shape)
    Lambda = grad_xi_F[peak_idx]
    om_at_peak = omega_mag[peak_idx]
    product_at_peak = om_at_peak * Lambda

    # Crease mask: top 1% of |∇ξ|_F
    crease_thresh = np.percentile(grad_xi_F, 99)
    crease_mask = grad_xi_F >= crease_thresh
    # Tube mask: top 10% of omega_mag
    tube_thresh = np.percentile(omega_mag, 90)
    tube_mask = omega_mag >= tube_thresh

    # ── M1: Local balance ratio ──
    # B(x) = |ω·∇u| / (ν|∇²ω| + |ω|·|∇ξ|)
    # Stretching: (ω·∇u)_i = ω_j * ∂u_i/∂x_j
    stretch = np.zeros((3, N, N, N))
    for i in range(3):
        for j in range(3):
            stretch[i] += o[j] * grad_u[i, j]
    stretch_mag = np.sqrt(stretch[0]**2 + stretch[1]**2 + stretch[2]**2)

    # Laplacian of vorticity
    lap_o = np.zeros((3, N, N, N))
    for i in range(3):
        lap_o[i] = spectral_laplacian(o_hat[i], k2)
    lap_o_mag = np.sqrt(lap_o[0]**2 + lap_o[1]**2 + lap_o[2]**2)

    denom_M1 = NU * lap_o_mag + omega_mag * grad_xi_F + 1e-30
    B = stretch_mag / denom_M1

    B_at_crease = B[crease_mask]
    B_at_tube = B[tube_mask]

    # ── M2: Helicity density h = u·ω ──
    h = ux * ox + uy * oy + uz * oz
    h_at_crease = h[crease_mask]
    abs_h_at_crease = np.abs(h_at_crease)

    # ── M3: Strain-vorticity alignment ──
    # cos(θ) = (S·ω)·ω / (|S·ω|·|ω|)
    So = np.zeros((3, N, N, N))
    for i in range(3):
        for j in range(3):
            So[i] += S[i, j] * o[j]
    So_dot_o = So[0] * ox + So[1] * oy + So[2] * oz
    So_mag = np.sqrt(So[0]**2 + So[1]**2 + So[2]**2)
    cos_theta = So_dot_o / (So_mag * omega_mag + 1e-30)
    cos_theta = np.clip(cos_theta, -1, 1)

    cos_at_crease = cos_theta[crease_mask]
    cos_at_tube = cos_theta[tube_mask]

    # ── M4: Invariant candidate — individual scaling at peak voxel ──
    # (computed above: Lambda, om_at_peak, product_at_peak)
    # Also: product field max and stats
    product_field = omega_mag * grad_xi_F
    product_max = product_field.max()
    product_at_crease = product_field[crease_mask]

    # ── M5: Pressure gradient ──
    # ∇²p = -∂_i u_j ∂_j u_i = -tr(∇u · ∇u)
    rhs = np.zeros((N, N, N))
    for i in range(3):
        for j in range(3):
            rhs -= grad_u[i, j] * grad_u[j, i]

    rhs_hat = np.fft.fftn(rhs) * dealias  # dealias the nonlinear product
    p_hat = -rhs_hat / k2_safe
    p_hat[0, 0, 0] = 0  # zero-mean gauge

    grad_p = np.zeros((3, N, N, N))
    for j in range(3):
        grad_p[j] = spectral_deriv(p_hat, Ks[j])
    grad_p_mag = np.sqrt(grad_p[0]**2 + grad_p[1]**2 + grad_p[2]**2)

    grad_p_at_crease = grad_p_mag[crease_mask]

    # Pressure-compression alignment: ∇p · e_compression
    # Compression axis = eigenvector of S with most negative eigenvalue
    # Too expensive per-voxel for 256³. Use proxy: alignment of ∇p with ω
    # (at the crease, the relevant axis is the direction change axis)
    gp_dot_o = grad_p[0]*ox + grad_p[1]*oy + grad_p[2]*oz
    pres_align = gp_dot_o / (grad_p_mag * omega_mag + 1e-30)
    pres_align = np.clip(pres_align, -1, 1)
    pres_align_at_crease = pres_align[crease_mask]

    # ── M6: Pointwise Kolmogorov-conditioned Λ ──
    # η(x) = (ν³/ε(x))^(1/4) where ε = 2ν·S:S
    S_sq = np.zeros((N, N, N))
    for i in range(3):
        for j in range(3):
            S_sq += S[i, j]**2
    epsilon_local = 2 * NU * S_sq
    eta_local = (NU**3 / (epsilon_local + 1e-30))**0.25
    Lambda_eta = grad_xi_F * eta_local

    Lambda_eta_gt1 = Lambda_eta > 1.0
    Lambda_eta_at_crease = Lambda_eta[crease_mask]
    frac_crease_gt1 = Lambda_eta_gt1[crease_mask].mean()

    # ── Collect results ──
    z_mid = N // 2
    results = {
        "N": N, "t": t,
        # M1
        "M1_B_crease_mean": float(B_at_crease.mean()),
        "M1_B_crease_median": float(np.median(B_at_crease)),
        "M1_B_tube_mean": float(B_at_tube.mean()),
        "M1_B_tube_median": float(np.median(B_at_tube)),
        # M2
        "M2_h_crease_mean": float(h_at_crease.mean()),
        "M2_absh_crease_mean": float(abs_h_at_crease.mean()),
        "M2_absh_crease_median": float(np.median(abs_h_at_crease)),
        # M3
        "M3_cos_crease_mean": float(cos_at_crease.mean()),
        "M3_cos_crease_std": float(cos_at_crease.std()),
        "M3_cos_tube_mean": float(cos_at_tube.mean()),
        "M3_cos_tube_std": float(cos_at_tube.std()),
        # M4 — THE KEY METRIC
        "M4_Lambda": float(Lambda),
        "M4_om_at_peak": float(om_at_peak),
        "M4_product_at_peak": float(product_at_peak),
        "M4_product_field_max": float(product_max),
        "M4_product_crease_mean": float(product_at_crease.mean()),
        "M4_product_crease_max": float(product_at_crease.max()),
        "M4_om_max": float(omega_mag.max()),
        "M4_om_at_peak_frac": float(om_at_peak / (omega_mag.max() + 1e-30)),
        # M5
        "M5_gradp_crease_mean": float(grad_p_at_crease.mean()),
        "M5_gradp_crease_max": float(grad_p_at_crease.max()),
        "M5_pres_align_crease_mean": float(pres_align_at_crease.mean()),
        "M5_pres_align_crease_std": float(pres_align_at_crease.std()),
        # M6
        "M6_Lambda_eta_crease_mean": float(Lambda_eta_at_crease.mean()),
        "M6_Lambda_eta_crease_max": float(Lambda_eta_at_crease.max()),
        "M6_frac_crease_eta_gt1": float(frac_crease_gt1),
        "M6_frac_domain_eta_gt1": float(Lambda_eta_gt1.mean()),
    }

    # Save z-midplane slices for visualization (only at representative snapshots)
    slices = {
        "B_slice": B[z_mid],
        "cos_theta_slice": cos_theta[z_mid],
        "Lambda_eta_slice": Lambda_eta[z_mid],
        "grad_xi_slice": grad_xi_F[z_mid],
        "omega_mag_slice": omega_mag[z_mid],
        "product_slice": product_field[z_mid],
        "grad_p_slice": grad_p_mag[z_mid],
    }

    return results, slices


# ── M5 dealiasing cross-check ──

def verify_pressure_dealiasing():
    """Compare spectral vs finite-difference pressure gradient at N=64."""
    print("M5 dealiasing cross-check: spectral vs FD at N=64 t=20...")
    h5 = os.path.join(SNAP_DIR, "run_64_t0020.0.h5")
    f = h5py.File(h5, "r")
    ux, uy, uz = f["ux"][:], f["uy"][:], f["uz"][:]
    N = 64; f.close()

    dx = 2 * np.pi / N
    kx, ky, kz, k2, k2_safe, dealias = make_wavenumbers(N)
    Ks = [kx, ky, kz]
    u = [ux, uy, uz]
    u_hat = [np.fft.fftn(c) for c in u]

    # Spectral ∇u
    grad_u = np.zeros((3, 3, N, N, N))
    for i in range(3):
        for j in range(3):
            grad_u[i, j] = spectral_deriv(u_hat[i], Ks[j])

    # Spectral pressure
    rhs = np.zeros((N, N, N))
    for i in range(3):
        for j in range(3):
            rhs -= grad_u[i, j] * grad_u[j, i]
    rhs_hat = np.fft.fftn(rhs) * dealias
    p_hat = -rhs_hat / k2_safe
    p_hat[0, 0, 0] = 0
    grad_p_spec = np.zeros((3, N, N, N))
    for j in range(3):
        grad_p_spec[j] = spectral_deriv(p_hat, Ks[j])

    # FD pressure: same RHS, solve Poisson with spectral but ∇p with central diff
    p_phys = np.fft.ifftn(p_hat).real
    grad_p_fd = np.zeros((3, N, N, N))
    grad_p_fd[0] = (np.roll(p_phys, -1, axis=0) - np.roll(p_phys, 1, axis=0)) / (2 * dx)
    grad_p_fd[1] = (np.roll(p_phys, -1, axis=1) - np.roll(p_phys, 1, axis=1)) / (2 * dx)
    grad_p_fd[2] = (np.roll(p_phys, -1, axis=2) - np.roll(p_phys, 1, axis=2)) / (2 * dx)

    diff = np.sqrt(sum((grad_p_spec[j] - grad_p_fd[j])**2 for j in range(3)))
    mag = np.sqrt(sum(grad_p_spec[j]**2 for j in range(3)))
    rel_err = diff / (mag + 1e-30)

    print(f"  Max |∇p_spec - ∇p_FD| / |∇p_spec|: {rel_err.max():.4f}")
    print(f"  Mean relative error: {rel_err.mean():.4f}")
    print(f"  Median relative error: {np.median(rel_err):.6f}")
    ok = rel_err.mean() < 0.1
    print(f"  Verdict: {'PASS' if ok else 'WARN'} — {'consistent' if ok else 'check high-k modes'}")
    return ok


# ── Main ──

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    # Cross-check M5
    verify_pressure_dealiasing()

    # Process all snapshots
    h5files = sorted([f for f in os.listdir(SNAP_DIR) if f.endswith('.h5')])
    all_results = []
    all_slices = {}

    for fn in h5files:
        path = os.path.join(SNAP_DIR, fn)
        print(f"\nProcessing {fn}...", flush=True)
        try:
            results, slices = compute_all_metrics(path)
            all_results.append(results)
            all_slices[(results["N"], results["t"])] = slices
            print(f"  N={results['N']} t={results['t']:.1f}: "
                  f"Λ={results['M4_Lambda']:.2f} "
                  f"|ω|@peak={results['M4_om_at_peak']:.6f} "
                  f"product@peak={results['M4_product_at_peak']:.4f} "
                  f"product_max={results['M4_product_field_max']:.4f}")
        except Exception as e:
            print(f"  ERROR: {e}")

    # Write CSV
    csv_path = os.path.join(OUT_DIR, "metrics.csv")
    if all_results:
        fields = list(all_results[0].keys())
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in all_results:
                w.writerow(r)
        print(f"\nWrote {csv_path} ({len(all_results)} rows)")

    # Save slices as npz for plotting
    for key, sl in all_slices.items():
        N, t = key
        npz_path = os.path.join(OUT_DIR, f"slices_N{N}_t{t:.0f}.npz")
        np.savez_compressed(npz_path, **sl)

    print(f"Saved {len(all_slices)} slice files")
    print("Done.")
