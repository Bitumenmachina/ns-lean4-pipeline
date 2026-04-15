#!/usr/bin/env python3
"""
V5 — Vortex Reconnection Direction-Field Test

Does the vorticity direction field develop a discontinuity during
reconnection faster than viscous dissipation can smooth it?

Pseudospectral FFT solver, RK4 adaptive, 2/3 dealiased.
Antiparallel vortex tubes with forced reconnection.

Kill metric K = τ_dir / τ_diss:
  K > 1 and growing → direction tearing outpaces dissipation
  K ≤ 1 → dissipation catches

Usage: python3 ns_v5_reconnection.py [--data-dir PATH]
       NS_DATA_DIR=/path/to/output python3 ns_v5_reconnection.py
"""

import argparse
import json
import os
import sqlite3
import sys
import time
from pathlib import Path

import numpy as np

BASE_DIR = Path(
    os.environ.get("NS_DATA_DIR", "./data/v5")
)
DB_PATH = BASE_DIR / "reconnection_v5.db"

# ============================================================================
# PSEUDOSPECTRAL 3D NS SOLVER
# ============================================================================

class PseudospectralNS3D:
    """
    3D incompressible NS on [0, 2π]³ periodic.
    Pseudospectral: derivatives exact in Fourier space.
    2/3 dealiasing. RK4 adaptive timestepping.
    """

    def __init__(self, N, nu):
        self.N = N
        self.nu = nu
        self.L = 2 * np.pi
        self.dx = self.L / N

        # Physical space grid
        x = np.linspace(0, self.L, N, endpoint=False)
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing='ij')

        # Wavenumber arrays for rfft (last axis halved)
        kx = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
        ky = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
        kz = np.fft.rfftfreq(N, d=1.0 / N).astype(np.float64)
        self.KX, self.KY, self.KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        self.K2 = self.KX**2 + self.KY**2 + self.KZ**2
        self.K2_safe = np.where(self.K2 == 0, 1.0, self.K2)

        # Dealiasing mask: 2/3 rule
        k_max = N // 3
        self.dealias = (
            (np.abs(self.KX) <= k_max) &
            (np.abs(self.KY) <= k_max) &
            (np.abs(self.KZ) <= k_max)
        ).astype(np.float64)

        # Velocity in Fourier space (complex, rfft shape)
        shape_k = (N, N, N // 2 + 1)
        self.ux_hat = np.zeros(shape_k, dtype=np.complex128)
        self.uy_hat = np.zeros(shape_k, dtype=np.complex128)
        self.uz_hat = np.zeros(shape_k, dtype=np.complex128)

        self.t = 0.0
        self.step_count = 0
        self._update_dt()

    def _update_dt(self):
        """Compute adaptive dt from CFL condition."""
        ux = np.fft.irfftn(self.ux_hat, axes=(0, 1, 2))
        uy = np.fft.irfftn(self.uy_hat, axes=(0, 1, 2))
        uz = np.fft.irfftn(self.uz_hat, axes=(0, 1, 2))
        u_max = max(np.max(np.abs(ux)), np.max(np.abs(uy)),
                    np.max(np.abs(uz)), 1e-10)
        self.dt = min(0.3 * self.dx / u_max, 0.01)

    def init_antiparallel_tubes(self):
        """Antiparallel vortex tubes with sinusoidal perturbation."""
        N = self.N
        d0 = np.pi / 2       # tube separation
        sigma0 = np.pi / 8   # core radius
        Gamma = 1.0           # circulation
        A_pert = d0 / 4       # perturbation amplitude

        cx = np.pi            # center of domain
        cy = np.pi

        # Vorticity field (only z-component for antiparallel tubes along z)
        omega_z = np.zeros((N, N, N), dtype=np.float64)

        for sign, offset_sign in [(+1, +1), (-1, -1)]:
            # Perturbation shifts tubes toward each other at z=π
            x_offset = offset_sign * (d0 / 2) - offset_sign * A_pert * np.cos(self.Z)

            r2 = (self.X - cx - x_offset)**2 + (self.Y - cy)**2
            omega_z += sign * Gamma / (np.pi * sigma0**2) * np.exp(-r2 / sigma0**2)

        # Get velocity from vorticity via Biot-Savart in Fourier space
        # û = (ik × ω̂) / k²
        # With ω = (0, 0, ω_z): ik × ω = (ik_y ω_z, -ik_x ω_z, 0) / k²
        omega_z_hat = np.fft.rfftn(omega_z, axes=(0, 1, 2))

        self.ux_hat = (1j * self.KY * omega_z_hat) / self.K2_safe
        self.uy_hat = (-1j * self.KX * omega_z_hat) / self.K2_safe
        self.uz_hat = np.zeros_like(self.ux_hat)

        # Zero mean
        self.ux_hat[0, 0, 0] = 0
        self.uy_hat[0, 0, 0] = 0

        # Ensure divergence-free
        self._project()
        self._update_dt()

    def _project(self):
        """Remove divergent component: û -= k(k·û)/k²"""
        kdotu = (self.KX * self.ux_hat + self.KY * self.uy_hat +
                 self.KZ * self.uz_hat)
        self.ux_hat -= self.KX * kdotu / self.K2_safe
        self.uy_hat -= self.KY * kdotu / self.K2_safe
        self.uz_hat -= self.KZ * kdotu / self.K2_safe

    def _rhs(self, ux_hat, uy_hat, uz_hat):
        """Compute dû/dt = -FFT(u·∇u) - νk²û, then project."""
        # To physical space
        ux = np.fft.irfftn(ux_hat, axes=(0, 1, 2))
        uy = np.fft.irfftn(uy_hat, axes=(0, 1, 2))
        uz = np.fft.irfftn(uz_hat, axes=(0, 1, 2))

        # Compute u_i * u_j products in physical space (6 unique)
        uxux = np.fft.rfftn(ux * ux, axes=(0, 1, 2)) * self.dealias
        uxuy = np.fft.rfftn(ux * uy, axes=(0, 1, 2)) * self.dealias
        uxuz = np.fft.rfftn(ux * uz, axes=(0, 1, 2)) * self.dealias
        uyuy = np.fft.rfftn(uy * uy, axes=(0, 1, 2)) * self.dealias
        uyuz = np.fft.rfftn(uy * uz, axes=(0, 1, 2)) * self.dealias
        uzuz = np.fft.rfftn(uz * uz, axes=(0, 1, 2)) * self.dealias

        # Convection: -∂_j(u_i u_j) in Fourier space
        conv_x = -(1j * self.KX * uxux + 1j * self.KY * uxuy + 1j * self.KZ * uxuz)
        conv_y = -(1j * self.KX * uxuy + 1j * self.KY * uyuy + 1j * self.KZ * uyuz)
        conv_z = -(1j * self.KX * uxuz + 1j * self.KY * uyuz + 1j * self.KZ * uzuz)

        # Viscous: -νk²û
        visc = -self.nu * self.K2

        dux = conv_x + visc * ux_hat
        duy = conv_y + visc * uy_hat
        duz = conv_z + visc * uz_hat

        # Pressure projection
        kdot = self.KX * dux + self.KY * duy + self.KZ * duz
        dux -= self.KX * kdot / self.K2_safe
        duy -= self.KY * kdot / self.K2_safe
        duz -= self.KZ * kdot / self.K2_safe

        return dux, duy, duz

    def step_rk4(self):
        """One RK4 step with adaptive CFL."""
        dt = self.dt
        u = (self.ux_hat.copy(), self.uy_hat.copy(), self.uz_hat.copy())

        # k1
        k1 = self._rhs(*u)
        # k2
        u2 = tuple(u[i] + 0.5 * dt * k1[i] for i in range(3))
        k2 = self._rhs(*u2)
        # k3
        u3 = tuple(u[i] + 0.5 * dt * k2[i] for i in range(3))
        k3 = self._rhs(*u3)
        # k4
        u4 = tuple(u[i] + dt * k3[i] for i in range(3))
        k4 = self._rhs(*u4)

        # Update
        for i, attr in enumerate(['ux_hat', 'uy_hat', 'uz_hat']):
            setattr(self, attr,
                    u[i] + (dt / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]))

        self.t += dt
        self.step_count += 1

        # Update dt every 10 steps
        if self.step_count % 10 == 0:
            self._update_dt()

    def get_physical(self):
        """Return velocity in physical space."""
        ux = np.fft.irfftn(self.ux_hat, axes=(0, 1, 2))
        uy = np.fft.irfftn(self.uy_hat, axes=(0, 1, 2))
        uz = np.fft.irfftn(self.uz_hat, axes=(0, 1, 2))
        return ux, uy, uz

    def compute_vorticity_spectral(self):
        """Compute vorticity in Fourier space: ω̂ = ik × û."""
        ox_hat = 1j * self.KY * self.uz_hat - 1j * self.KZ * self.uy_hat
        oy_hat = 1j * self.KZ * self.ux_hat - 1j * self.KX * self.uz_hat
        oz_hat = 1j * self.KX * self.uy_hat - 1j * self.KY * self.ux_hat
        return ox_hat, oy_hat, oz_hat


# ============================================================================
# PSEUDOSPECTRAL 2D NS SOLVER (control)
# ============================================================================

class PseudospectralNS2D:
    """2D incompressible NS on [0, 2π]² periodic. Pseudospectral."""

    def __init__(self, N, nu):
        self.N = N
        self.nu = nu
        self.L = 2 * np.pi
        self.dx = self.L / N

        x = np.linspace(0, self.L, N, endpoint=False)
        self.X, self.Y = np.meshgrid(x, x, indexing='ij')

        kx = np.fft.fftfreq(N, d=1.0 / N).astype(np.float64)
        ky = np.fft.rfftfreq(N, d=1.0 / N).astype(np.float64)
        self.KX, self.KY = np.meshgrid(kx, ky, indexing='ij')
        self.K2 = self.KX**2 + self.KY**2
        self.K2_safe = np.where(self.K2 == 0, 1.0, self.K2)

        k_max = N // 3
        self.dealias = (
            (np.abs(self.KX) <= k_max) & (np.abs(self.KY) <= k_max)
        ).astype(np.float64)

        shape_k = (N, N // 2 + 1)
        self.ux_hat = np.zeros(shape_k, dtype=np.complex128)
        self.uy_hat = np.zeros(shape_k, dtype=np.complex128)
        self.t = 0.0
        self.step_count = 0
        self._update_dt()

    def _update_dt(self):
        ux = np.fft.irfft2(self.ux_hat)
        uy = np.fft.irfft2(self.uy_hat)
        u_max = max(np.max(np.abs(ux)), np.max(np.abs(uy)), 1e-10)
        self.dt = min(0.3 * self.dx / u_max, 0.01)

    def init_antiparallel_vortices(self):
        """2D antiparallel vortex pair."""
        d0 = np.pi / 2
        sigma0 = np.pi / 8
        Gamma = 1.0
        cx, cy = np.pi, np.pi

        omega = np.zeros((self.N, self.N), dtype=np.float64)
        for sign, xoff in [(+1, +d0/2), (-1, -d0/2)]:
            r2 = (self.X - cx - xoff)**2 + (self.Y - cy)**2
            omega += sign * Gamma / (np.pi * sigma0**2) * np.exp(-r2 / sigma0**2)

        # Stream function: ψ̂ = -ω̂/k²
        omega_hat = np.fft.rfft2(omega)
        psi_hat = -omega_hat / self.K2_safe
        psi_hat[0, 0] = 0

        # u = (∂ψ/∂y, -∂ψ/∂x)
        self.ux_hat = 1j * self.KY * psi_hat
        self.uy_hat = -1j * self.KX * psi_hat
        self._update_dt()

    def _rhs(self, ux_hat, uy_hat):
        ux = np.fft.irfft2(ux_hat)
        uy = np.fft.irfft2(uy_hat)

        uxux = np.fft.rfft2(ux * ux) * self.dealias
        uxuy = np.fft.rfft2(ux * uy) * self.dealias
        uyuy = np.fft.rfft2(uy * uy) * self.dealias

        conv_x = -(1j * self.KX * uxux + 1j * self.KY * uxuy)
        conv_y = -(1j * self.KX * uxuy + 1j * self.KY * uyuy)

        visc = -self.nu * self.K2
        dux = conv_x + visc * ux_hat
        duy = conv_y + visc * uy_hat

        kdot = self.KX * dux + self.KY * duy
        dux -= self.KX * kdot / self.K2_safe
        duy -= self.KY * kdot / self.K2_safe
        return dux, duy

    def step_rk4(self):
        dt = self.dt
        u = (self.ux_hat.copy(), self.uy_hat.copy())
        k1 = self._rhs(*u)
        k2 = self._rhs(u[0] + 0.5*dt*k1[0], u[1] + 0.5*dt*k1[1])
        k3 = self._rhs(u[0] + 0.5*dt*k2[0], u[1] + 0.5*dt*k2[1])
        k4 = self._rhs(u[0] + dt*k3[0], u[1] + dt*k3[1])
        self.ux_hat = u[0] + (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        self.uy_hat = u[1] + (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        self.t += dt
        self.step_count += 1
        if self.step_count % 10 == 0:
            self._update_dt()

    def compute_vorticity_spectral(self):
        """2D vorticity: scalar ω = ∂uy/∂x - ∂ux/∂y"""
        return 1j * self.KX * self.uy_hat - 1j * self.KY * self.ux_hat


# ============================================================================
# DIAGNOSTICS: DIRECTION FIELD AND KILL METRIC
# ============================================================================

def compute_diagnostics_3d(solver):
    """Compute all V5 observables from the 3D solver state."""
    N = solver.N

    # Vorticity in Fourier space
    ox_hat, oy_hat, oz_hat = solver.compute_vorticity_spectral()

    # Physical space vorticity
    ox = np.fft.irfftn(ox_hat, axes=(0, 1, 2))
    oy = np.fft.irfftn(oy_hat, axes=(0, 1, 2))
    oz = np.fft.irfftn(oz_hat, axes=(0, 1, 2))

    omega_mag = np.sqrt(ox**2 + oy**2 + oz**2)
    omega_max = float(np.max(omega_mag))

    # Enstrophy
    enstrophy = float(np.mean(omega_mag**2)) * (2 * np.pi)**3

    # Velocity for strain rate
    ux, uy, uz = solver.get_physical()

    # Strain rate tensor S_ij = (∂u_i/∂x_j + ∂u_j/∂x_i) / 2
    # Dissipation: 2ν ∫ |S|² dV = 2ν ∫ S_ij S_ij dV
    # In Fourier: ∂u_i/∂x_j → ik_j û_i
    S2_sum = np.zeros((N, N, N), dtype=np.float64)
    u_hats = [solver.ux_hat, solver.uy_hat, solver.uz_hat]
    Ks = [solver.KX, solver.KY, solver.KZ]
    for i in range(3):
        for j in range(3):
            du_i_dx_j = np.fft.irfftn(1j * Ks[j] * u_hats[i], axes=(0, 1, 2))
            du_j_dx_i = np.fft.irfftn(1j * Ks[i] * u_hats[j], axes=(0, 1, 2))
            S_ij = 0.5 * (du_i_dx_j + du_j_dx_i)
            S2_sum += S_ij**2

    dissipation = float(2 * solver.nu * np.mean(S2_sum) * (2 * np.pi)**3)

    # Velocity gradient norm |A|
    A2_sum = np.zeros((N, N, N), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            du_dx = np.fft.irfftn(1j * Ks[j] * u_hats[i], axes=(0, 1, 2))
            A2_sum += du_dx**2
    A_max = float(np.sqrt(np.max(A2_sum)))

    # ---- DIRECTION FIELD ξ = ω/|ω| ----
    # Threshold: only where |ω| > 1% of max
    eps_omega = 0.01 * omega_max if omega_max > 1e-10 else 1e-10
    valid = omega_mag > eps_omega
    valid_count = int(np.sum(valid))

    if valid_count < 100:
        # Not enough vorticity to compute direction field
        return {
            "Lambda_max": 0.0, "enstrophy": enstrophy,
            "dissipation": dissipation, "R_local": omega_max,
            "omega_max": omega_max, "A_max": A_max,
            "valid_count": valid_count,
        }

    # Direction field ξ_i = ω_i / |ω|
    safe_mag = np.maximum(omega_mag, eps_omega)
    xi_x = ox / safe_mag
    xi_y = oy / safe_mag
    xi_z = oz / safe_mag

    # ∇ξ: gradient of direction field
    # ∂ξᵢ/∂xⱼ = (∂ωᵢ/∂xⱼ)/|ω| - ξᵢ(ξₖ ∂ωₖ/∂xⱼ)/|ω|
    # All derivatives spectral
    o_hats = [ox_hat, oy_hat, oz_hat]
    xis = [xi_x, xi_y, xi_z]

    grad_xi_F2 = np.zeros((N, N, N), dtype=np.float64)  # Frobenius norm squared

    for j in range(3):
        # ∂ω_k/∂x_j for all k
        d_omega_dx_j = []
        for k in range(3):
            d_omega_dx_j.append(np.fft.irfftn(1j * Ks[j] * o_hats[k], axes=(0, 1, 2)))

        # ξ_k * ∂ω_k/∂x_j (sum over k)
        xi_dot_grad_omega_j = sum(xis[k] * d_omega_dx_j[k] for k in range(3))

        for i in range(3):
            # ∂ξᵢ/∂xⱼ = (∂ωᵢ/∂xⱼ - ξᵢ * ξₖ∂ωₖ/∂xⱼ) / |ω|
            dxi_i_dxj = (d_omega_dx_j[i] - xis[i] * xi_dot_grad_omega_j) / safe_mag
            grad_xi_F2 += dxi_i_dxj**2

    grad_xi_F = np.sqrt(grad_xi_F2)

    # Λ = max |∇ξ|_F where |ω| > threshold
    Lambda_max = float(np.max(grad_xi_F[valid]))

    return {
        "Lambda_max": Lambda_max,
        "enstrophy": enstrophy,
        "dissipation": dissipation,
        "R_local": omega_max,
        "omega_max": omega_max,
        "A_max": A_max,
        "valid_count": valid_count,
    }


def compute_diagnostics_2d(solver):
    """Compute 2D diagnostics. No stretching, scalar vorticity."""
    omega_hat = solver.compute_vorticity_spectral()
    omega = np.fft.irfft2(omega_hat)
    omega_mag = np.abs(omega)
    omega_max = float(np.max(omega_mag))

    enstrophy = float(np.mean(omega**2)) * (2 * np.pi)**2

    # Dissipation
    ux = np.fft.irfft2(solver.ux_hat)
    uy = np.fft.irfft2(solver.uy_hat)
    N = solver.N

    dux_dx = np.fft.irfft2(1j * solver.KX * solver.ux_hat)
    dux_dy = np.fft.irfft2(1j * solver.KY * solver.ux_hat)
    duy_dx = np.fft.irfft2(1j * solver.KX * solver.uy_hat)
    duy_dy = np.fft.irfft2(1j * solver.KY * solver.uy_hat)

    S2 = (0.5*(dux_dx+dux_dx))**2 + (0.5*(dux_dy+duy_dx))**2 + \
         (0.5*(duy_dx+dux_dy))**2 + (0.5*(duy_dy+duy_dy))**2
    dissipation = float(2 * solver.nu * np.mean(S2) * (2 * np.pi)**2)

    A_max = float(np.sqrt(np.max(dux_dx**2 + dux_dy**2 + duy_dx**2 + duy_dy**2)))

    # 2D "direction field": sign of vorticity (scalar, ±1)
    # Λ in 2D: max |∇(sign(ω))| — but sign is discontinuous at ω=0
    # Use continuous proxy: ξ = ω/|ω|, gradient = |∇ω|/|ω| - ...
    # For 2D scalar: ∂(ω/|ω|)/∂x = ∂ω/∂x / |ω| * (1 - (ω/|ω|)²) = 0 when |ω|>0
    # Since ω/|ω| = ±1 (constant where defined), ∇ξ = 0 everywhere except ω=0
    # So Lambda = 0 in 2D (the direction field is trivially constant)
    # Instead, compute gradient of vorticity magnitude as proxy
    eps_omega = 0.01 * omega_max if omega_max > 1e-10 else 1e-10
    valid = omega_mag > eps_omega

    # Gradient of omega
    domega_dx = np.fft.irfft2(1j * solver.KX * omega_hat)
    domega_dy = np.fft.irfft2(1j * solver.KY * omega_hat)
    grad_omega = np.sqrt(domega_dx**2 + domega_dy**2)

    # Normalized gradient: |∇ω|/|ω| (measures how fast direction could change)
    Lambda_proxy = np.where(valid, grad_omega / omega_mag, 0.0)
    Lambda_max = float(np.max(Lambda_proxy[valid])) if np.any(valid) else 0.0

    return {
        "Lambda_max": Lambda_max,
        "enstrophy": enstrophy,
        "dissipation": dissipation,
        "R_local": omega_max,
        "omega_max": omega_max,
        "A_max": A_max,
        "valid_count": int(np.sum(valid)),
    }


# ============================================================================
# DATABASE
# ============================================================================

def init_db():
    conn = sqlite3.connect(str(DB_PATH))
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS runs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        dimension INTEGER, Re REAL, N INTEGER, nu REAL,
        total_steps INTEGER, elapsed_seconds REAL,
        label TEXT, timestamp TEXT
    )""")
    c.execute("""CREATE TABLE IF NOT EXISTS timeseries (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        run_id INTEGER, step INTEGER, time_phys REAL,
        Lambda_max REAL, enstrophy REAL, dissipation REAL,
        R_local REAL, omega_max REAL, A_max REAL,
        tau_dir REAL, tau_diss REAL, K_metric REAL,
        valid_count INTEGER,
        FOREIGN KEY (run_id) REFERENCES runs(id)
    )""")
    c.execute("""CREATE TABLE IF NOT EXISTS test_results (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        test_name TEXT, test_number INTEGER,
        passed INTEGER, value REAL, details TEXT
    )""")
    conn.commit()
    return conn


# ============================================================================
# RUN EXECUTION
# ============================================================================

def run_experiment(conn, dimension, Re, N, label, max_time=20.0, diag_interval=0.1):
    """Run one simulation."""
    nu = (2 * np.pi) / Re
    print(f"\n  {label}: {dimension}D Re={Re} N={N} nu={nu:.6f}")

    if dimension == 3:
        solver = PseudospectralNS3D(N, nu)
        solver.init_antiparallel_tubes()
    else:
        solver = PseudospectralNS2D(N, nu)
        solver.init_antiparallel_vortices()

    c = conn.cursor()
    c.execute("""INSERT INTO runs (dimension, Re, N, nu, total_steps, label, timestamp)
                 VALUES (?, ?, ?, ?, 0, ?, datetime('now'))""",
              (dimension, Re, N, nu, label))
    run_id = c.lastrowid

    ts_rows = []
    prev_Lambda = 0.0
    prev_t = 0.0
    t_start = time.time()
    next_diag = 0.0
    total_steps = 0

    while solver.t < max_time:
        solver.step_rk4()
        total_steps += 1

        if solver.t >= next_diag:
            next_diag = solver.t + diag_interval

            if dimension == 3:
                diag = compute_diagnostics_3d(solver)
            else:
                diag = compute_diagnostics_2d(solver)

            # Kill metric
            dt_diag = solver.t - prev_t if prev_t > 0 else diag_interval
            tau_dir = (diag["Lambda_max"] - prev_Lambda) / max(dt_diag, 1e-10)
            tau_diss = diag["dissipation"] / max(diag["enstrophy"], 1e-15)
            K = abs(tau_dir) / max(abs(tau_diss), 1e-15)

            ts_rows.append((
                run_id, total_steps, solver.t,
                diag["Lambda_max"], diag["enstrophy"], diag["dissipation"],
                diag["R_local"], diag["omega_max"], diag["A_max"],
                tau_dir, tau_diss, K, diag["valid_count"],
            ))

            prev_Lambda = diag["Lambda_max"]
            prev_t = solver.t

            # Incremental commit every 50 diagnostics (~5 sim time units)
            if len(ts_rows) % 50 == 0 and ts_rows:
                c.executemany("""INSERT INTO timeseries
                    (run_id, step, time_phys, Lambda_max, enstrophy, dissipation,
                     R_local, omega_max, A_max, tau_dir, tau_diss, K_metric, valid_count)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", ts_rows)
                conn.commit()
                flushed = len(ts_rows)
                ts_rows.clear()
                print(f"    [checkpoint] {flushed} rows committed at t={solver.t:.2f}")

            # Progress every 2.0 time units
            if len(ts_rows) % max(1, int(2.0 / diag_interval)) == 0:
                print(f"    t={solver.t:.2f}/{max_time:.0f}  Λ={diag['Lambda_max']:.2f}  "
                      f"Ω={diag['enstrophy']:.4f}  ε={diag['dissipation']:.4f}  "
                      f"|ω|={diag['omega_max']:.2f}  K={K:.4f}  "
                      f"dt={solver.dt:.5f}")

            # Blowup check
            if diag["A_max"] > 1e8 or np.isnan(diag["A_max"]):
                print(f"    BLOWUP at t={solver.t:.4f}")
                break

    elapsed = time.time() - t_start

    if ts_rows:
        c.executemany("""INSERT INTO timeseries
            (run_id, step, time_phys, Lambda_max, enstrophy, dissipation,
             R_local, omega_max, A_max, tau_dir, tau_diss, K_metric, valid_count)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", ts_rows)

    c.execute("UPDATE runs SET total_steps = ?, elapsed_seconds = ? WHERE id = ?",
              (total_steps, elapsed, run_id))
    conn.commit()

    print(f"    Done: {elapsed:.1f}s, {total_steps} steps, {len(ts_rows)} diagnostics")
    return run_id


# ============================================================================
# TEST EVALUATION
# ============================================================================

def evaluate_tests(conn):
    c = conn.cursor()
    print("\n" + "=" * 60)
    print("TEST EVALUATION")
    print("=" * 60)

    # Get all runs
    c.execute("SELECT id, dimension, Re, N FROM runs ORDER BY dimension, Re, N")
    all_runs = c.fetchall()
    runs_3d = [(r[0], r[2], r[3]) for r in all_runs if r[1] == 3]
    runs_2d = [(r[0], r[2], r[3]) for r in all_runs if r[1] == 2]

    # ---- T1: Λ increases during reconnection ----
    details = ""
    peak_lambdas = {}
    for run_id, Re, N in runs_3d:
        c.execute("""SELECT time_phys, Lambda_max, enstrophy FROM timeseries
                     WHERE run_id = ? ORDER BY step""", (run_id,))
        rows = c.fetchall()
        if not rows:
            continue
        L_arr = [r[1] for r in rows]
        peak_L = max(L_arr)
        peak_idx = L_arr.index(peak_L)
        peak_t = rows[peak_idx][0]
        Z_at_peak = rows[peak_idx][2]
        peak_lambdas[(Re, N)] = peak_L
        details += f"  Re={Re} N={N}: Λ_peak={peak_L:.2f} at t={peak_t:.2f}, Z={Z_at_peak:.4f}\n"

    t1_pass = any(v > 10.0 for v in peak_lambdas.values())
    c.execute("INSERT INTO test_results VALUES (NULL, 'T1: Λ peaks during reconnection', 1, ?, ?, ?)",
              (1 if t1_pass else 0, max(peak_lambdas.values()) if peak_lambdas else 0, details))

    # ---- T2: Λ scales with Re ----
    details = ""
    n128_data = [(Re, peak_lambdas.get((Re, 128), 0)) for Re, N in peak_lambdas if N == 128]
    n128_data = [(Re, L) for Re, L in sorted(set(n128_data)) if L > 0]
    if len(n128_data) >= 3:
        re_arr = np.array([d[0] for d in n128_data])
        L_arr = np.array([d[1] for d in n128_data])
        log_re = np.log10(re_arr)
        log_L = np.log10(np.clip(L_arr, 1e-10, None))
        alpha = np.polyfit(log_re, log_L, 1)[0]
        details = f"  Fit: Λ_peak ~ Re^{alpha:.3f}\n"
        for Re, L in n128_data:
            details += f"  Re={Re}: Λ_peak={L:.2f}\n"
    else:
        alpha = 0

    t2_pass = alpha > 0.5
    c.execute("INSERT INTO test_results VALUES (NULL, 'T2: Λ scales with Re (α>0.5)', 2, ?, ?, ?)",
              (1 if t2_pass else 0, alpha, details))

    # ---- T3: Kill metric K > 1 ----
    details = ""
    max_K_by_Re = {}
    for run_id, Re, N in runs_3d:
        c.execute("SELECT time_phys, K_metric FROM timeseries WHERE run_id = ? ORDER BY step",
                  (run_id,))
        rows = c.fetchall()
        if rows:
            K_vals = [r[1] for r in rows if r[1] is not None]
            if K_vals:
                peak_K = max(K_vals)
                peak_t = rows[K_vals.index(peak_K)][0]
                max_K_by_Re[(Re, N)] = peak_K
                details += f"  Re={Re} N={N}: K_peak={peak_K:.4f} at t={peak_t:.2f}\n"

    t3_pass = any(v > 1.0 for v in max_K_by_Re.values())
    c.execute("INSERT INTO test_results VALUES (NULL, 'T3: Kill metric K>1', 3, ?, ?, ?)",
              (1 if t3_pass else 0, max(max_K_by_Re.values()) if max_K_by_Re else 0, details))

    # ---- T4: Resolution convergence ----
    details = ""
    converged = True
    for Re_target in [1600, 3200]:
        N_data = []
        for Re, N in peak_lambdas:
            if abs(Re - Re_target) < 1:
                N_data.append((N, peak_lambdas[(Re, N)]))
        N_data.sort()
        if len(N_data) >= 2:
            details += f"  Re={Re_target}: "
            for N, L in N_data:
                details += f"N={N}→Λ={L:.2f}  "
            details += "\n"
            # Check: does Λ diverge or converge?
            L_vals = [d[1] for d in N_data]
            if len(L_vals) >= 2:
                ratio = L_vals[-1] / max(L_vals[-2], 1e-10)
                if ratio > 3.0:
                    details += f"    Ratio={ratio:.2f} — DIVERGING (potential singularity)\n"
                elif ratio > 1.5:
                    details += f"    Ratio={ratio:.2f} — not yet converged\n"
                    converged = False
                else:
                    details += f"    Ratio={ratio:.2f} — converging\n"

    c.execute("INSERT INTO test_results VALUES (NULL, 'T4: Λ converges with resolution', 4, ?, ?, ?)",
              (1 if converged else 0, 0, details))

    # ---- T5: 2D control ----
    details = ""
    max_K_2d = 0
    for run_id, Re, N in runs_2d:
        c.execute("SELECT time_phys, K_metric, Lambda_max FROM timeseries WHERE run_id = ? ORDER BY step",
                  (run_id,))
        rows = c.fetchall()
        if rows:
            K_vals = [r[1] for r in rows if r[1] is not None]
            L_vals = [r[2] for r in rows if r[2] is not None]
            if K_vals:
                pk = max(K_vals)
                max_K_2d = max(max_K_2d, pk)
                details += f"  2D Re={Re}: K_peak={pk:.4f}, Λ_peak={max(L_vals):.2f}\n"

    t5_pass = max_K_2d < 1.0
    c.execute("INSERT INTO test_results VALUES (NULL, 'T5: 2D control (K<1)', 5, ?, ?, ?)",
              (1 if t5_pass else 0, max_K_2d, details))

    conn.commit()

    # Print
    c.execute("SELECT test_name, test_number, passed, value, details FROM test_results ORDER BY test_number")
    for name, num, p, val, det in c.fetchall():
        status = "PASS" if p else "FAIL"
        print(f"\n  {name} — {status} (value={val:.4f})")
        for line in (det or "").split('\n')[:6]:
            if line.strip():
                print(f"    {line.strip()}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 60)
    print("V5 — VORTEX RECONNECTION DIRECTION-FIELD TEST")
    print("Pseudospectral FFT, RK4, Antiparallel Tubes")
    print("=" * 60)

    BASE_DIR.mkdir(parents=True, exist_ok=True)
    if DB_PATH.exists():
        DB_PATH.unlink()
    conn = init_db()

    total_start = time.time()

    # ===== 2D CONTROL (runs first) =====
    print("\n" + "=" * 60)
    print("2D CONTROL")
    print("=" * 60)
    for Re in [400, 1600]:
        run_experiment(conn, 2, Re, 128, f"2D_Re{Re}", max_time=15.0, diag_interval=0.1)

    # ===== 3D RE SWEEP at N=128 =====
    print("\n" + "=" * 60)
    print("3D RE SWEEP (N=128)")
    print("=" * 60)
    for Re in [400, 800, 1600, 3200, 6400]:
        run_experiment(conn, 3, Re, 128, f"3D_Re{Re}_N128", max_time=20.0, diag_interval=0.1)

    # ===== RESOLUTION CONVERGENCE at N=64 =====
    print("\n" + "=" * 60)
    print("RESOLUTION CONVERGENCE (N=64)")
    print("=" * 60)
    for Re in [1600, 3200]:
        run_experiment(conn, 3, Re, 64, f"3D_Re{Re}_N64", max_time=20.0, diag_interval=0.1)

    # ===== EVALUATE =====
    evaluate_tests(conn)

    total_elapsed = time.time() - total_start
    c = conn.cursor()
    c.execute("SELECT COUNT(*) FROM test_results WHERE passed = 1")
    passed = c.fetchone()[0]
    c.execute("SELECT COUNT(*) FROM test_results")
    total = c.fetchone()[0]

    print("\n" + "=" * 60)
    print(f"V5 COMPLETE: {passed}/{total} tests, {total_elapsed:.0f}s ({total_elapsed/3600:.1f} hours)")
    print(f"Database: {DB_PATH}")
    print("=" * 60)

    conn.close()
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="V5 Vortex Reconnection Direction-Field Test")
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="Output directory (default: $NS_DATA_DIR or ./data/v5)")
    args = parser.parse_args()
    if args.data_dir is not None:
        BASE_DIR = args.data_dir
        DB_PATH = BASE_DIR / "reconnection_v5.db"
    sys.exit(main())
