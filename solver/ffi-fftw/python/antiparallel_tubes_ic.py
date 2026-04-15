"""Antiparallel Gaussian vortex tubes IC (Roediger-Moriarty 2026 §4.1).

Domain: [0, 2π]³ periodic
Two antiparallel vortex tubes aligned along x-axis, separated in y,
with sinusoidal z-perturbation.

Parameters:
    Γ = 1          circulation
    σ = 0.3        core radius
    d = π/2        separation in y-direction
    A = 0.1        sinusoidal perturbation amplitude in z
"""
import numpy as np


def antiparallel_tubes_ic(n):
    """Generate velocity field for antiparallel vortex tubes.

    Each tube is a Gaussian vortex filament:
        ω_x = ±(Γ / (π σ²)) exp(-r²/σ²)

    where r is distance from filament centerline. The velocity field
    is obtained by inverting curl(ω) = -∇²u in Fourier space.

    Returns (ux, uy, uz) as (n, n, n) float64 arrays.
    """
    gamma = 1.0
    sigma = 0.3
    d = np.pi / 2.0
    amp = 0.1

    dx = 2 * np.pi / n
    x = np.arange(n) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    # Tube centerlines at y = π ± d/2, with z-perturbation
    y_center = np.pi
    y1 = y_center + d / 2.0
    y2 = y_center - d / 2.0

    # Sinusoidal perturbation of z-position along x
    z_center = np.pi
    z1 = z_center + amp * np.sin(X)
    z2 = z_center - amp * np.sin(X)

    # Periodic distance in y and z
    def pdist(a, b):
        """Periodic distance on [0, 2π]."""
        d = a - b
        return d - 2 * np.pi * np.round(d / (2 * np.pi))

    dy1 = pdist(Y, y1)
    dz1 = pdist(Z, z1)
    r1_sq = dy1**2 + dz1**2

    dy2 = pdist(Y, y2)
    dz2 = pdist(Z, z2)
    r2_sq = dy2**2 + dz2**2

    # Vorticity: omega_x component only (tubes along x-axis)
    # Tube 1: +Γ, Tube 2: -Γ (antiparallel)
    strength = gamma / (np.pi * sigma**2)
    omega_x = strength * (np.exp(-r1_sq / sigma**2) - np.exp(-r2_sq / sigma**2))
    omega_y = np.zeros_like(X)
    omega_z = np.zeros_like(X)

    # Invert curl: given ω, find u such that ∇×u = ω, ∇·u = 0
    # In Fourier space: û = (ik × ω̂) / k²
    kx = np.fft.fftfreq(n, d=1.0 / n)
    ky = np.fft.fftfreq(n, d=1.0 / n)
    kz = np.fft.rfftfreq(n, d=1.0 / n)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    K2 = KX**2 + KY**2 + KZ**2
    K2[0, 0, 0] = 1.0  # avoid division by zero

    ox_hat = np.fft.rfftn(omega_x)
    oy_hat = np.fft.rfftn(omega_y)
    oz_hat = np.fft.rfftn(omega_z)

    # û = (ik × ω̂) / k²
    # But we want u from ω = ∇×u, so û = -(ik × ω̂) / k²
    # Actually: ∇×u = ω → ik×û = ω̂ → û = (ik × ω̂) / (-k²)
    # Wait, let's be careful:
    # ∇×u = ω means i k × û = ω̂
    # So û = (ω̂ × ik) / k²  ... no.
    # ik × û = ω̂
    # Taking cross product with ik: ik × (ik × û) = ik × ω̂
    # ik × (ik × û) = ik(ik·û) - û(ik·ik) = -k² û  (since ∇·u = 0 → ik·û = 0)
    # So -k² û = ik × ω̂
    # û = -(ik × ω̂) / k²
    i = 1j
    cross_x = i * KY * oz_hat - i * KZ * oy_hat
    cross_y = i * KZ * ox_hat - i * KX * oz_hat
    cross_z = i * KX * oy_hat - i * KY * ox_hat

    ux_hat = -cross_x / K2
    uy_hat = -cross_y / K2
    uz_hat = -cross_z / K2

    # Zero mean flow
    ux_hat[0, 0, 0] = 0.0
    uy_hat[0, 0, 0] = 0.0
    uz_hat[0, 0, 0] = 0.0

    ux = np.fft.irfftn(ux_hat, s=(n, n, n), axes=(0, 1, 2))
    uy = np.fft.irfftn(uy_hat, s=(n, n, n), axes=(0, 1, 2))
    uz = np.fft.irfftn(uz_hat, s=(n, n, n), axes=(0, 1, 2))

    return ux, uy, uz


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    ux, uy, uz = antiparallel_tubes_ic(n)
    E = 0.5 * np.mean(ux**2 + uy**2 + uz**2)
    Om = np.mean(ux**2 + uy**2 + uz**2)  # rough proxy
    print(f"N={n}: E={E:.6f}, max|u|={np.sqrt(ux**2+uy**2+uz**2).max():.4f}")
    print(f"  ux range: [{ux.min():.4f}, {ux.max():.4f}]")
    print(f"  uy range: [{uy.min():.4f}, {uy.max():.4f}]")
    print(f"  uz range: [{uz.min():.4f}, {uz.max():.4f}]")
