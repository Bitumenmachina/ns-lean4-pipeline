//! Minimal pseudospectral Navier-Stokes kernel for 3D periodic domains.
//!
//! Rotational form: du/dt = P(ω × u) - ν k² u
//! Time integration: Integrating Factor RK4 (IFRK4)
//! Dealiasing: 2/3 rule
//! Domain: [0, 2π]³ periodic

use num_complex::Complex64;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use numpy::{PyReadonlyArray3, PyUntypedArrayMethods, PyArrayMethods, IntoPyArray};

/// Free-standing forward r2c FFT using pre-existing plan.
/// `phys_scratch` and `spec_scratch` are the working buffers for FFTW.
unsafe fn do_r2c(
    plan: fftw_sys::fftw_plan,
    src: &[f64],
    dst: &mut [Complex64],
    phys_scratch: &mut [f64],
    spec_scratch: &mut [Complex64],
) {
    phys_scratch[..src.len()].copy_from_slice(src);
    fftw_sys::fftw_execute_dft_r2c(
        plan,
        phys_scratch.as_mut_ptr(),
        spec_scratch.as_mut_ptr(),
    );
    dst.copy_from_slice(&spec_scratch[..dst.len()]);
}

/// Free-standing inverse c2r FFT using pre-existing plan. Normalizes by 1/N³.
/// c2r destroys input, so we copy into spec_scratch first.
unsafe fn do_c2r(
    plan: fftw_sys::fftw_plan,
    src: &[Complex64],
    dst: &mut [f64],
    phys_scratch: &mut [f64],
    spec_scratch: &mut [Complex64],
    norm: f64,
) {
    spec_scratch[..src.len()].copy_from_slice(src);
    fftw_sys::fftw_execute_dft_c2r(
        plan,
        spec_scratch.as_mut_ptr(),
        phys_scratch.as_mut_ptr(),
    );
    for i in 0..dst.len() {
        dst[i] = phys_scratch[i] * norm;
    }
}

/// Pseudospectral NS solver state and FFTW plans.
///
/// unsendable: fftw_plan is a raw pointer, not thread-safe.
/// This is fine — Python holds the GIL during all calls.
#[pyclass(unsendable)]
pub struct NSKernel {
    n: usize,
    nzc: usize,
    spec_len: usize,
    phys_len: usize,
    nu: f64,
    dt: f64,
    t: f64,
    step_count: usize,

    kx: Vec<f64>,
    ky: Vec<f64>,
    kz: Vec<f64>,

    k2: Vec<f64>,
    k2_safe: Vec<f64>,
    dealias: Vec<f64>,
    e_half: Vec<f64>,
    e_full: Vec<f64>,

    // Velocity state in Fourier space
    ux_hat: Vec<Complex64>,
    uy_hat: Vec<Complex64>,
    uz_hat: Vec<Complex64>,

    // FFTW plans
    plan_r2c: fftw_sys::fftw_plan,
    plan_c2r: fftw_sys::fftw_plan,

    // Planning buffers (kept alive for alignment consistency)
    _plan_phys_buf: Vec<f64>,
    _plan_spec_buf: Vec<Complex64>,

    // Scratch buffers for FFT I/O (shared across all FFT calls)
    fft_phys: Vec<f64>,
    fft_spec: Vec<Complex64>,

    // 6 physical-space scratch arrays for nonlinear term
    phys_scratch: [Vec<f64>; 6],

    // 3 spectral scratch for nonlinear RHS output
    nl_x: Vec<Complex64>,
    nl_y: Vec<Complex64>,
    nl_z: Vec<Complex64>,

    // 3 spectral scratch for vorticity (avoids clone in nonlinear_rhs)
    omega_x: Vec<Complex64>,
    omega_y: Vec<Complex64>,
    omega_z: Vec<Complex64>,

    // RK4 scratch
    u0_x: Vec<Complex64>,
    u0_y: Vec<Complex64>,
    u0_z: Vec<Complex64>,
    acc_x: Vec<Complex64>,
    acc_y: Vec<Complex64>,
    acc_z: Vec<Complex64>,
}

impl Drop for NSKernel {
    fn drop(&mut self) {
        unsafe {
            if !self.plan_r2c.is_null() {
                fftw_sys::fftw_destroy_plan(self.plan_r2c);
            }
            if !self.plan_c2r.is_null() {
                fftw_sys::fftw_destroy_plan(self.plan_c2r);
            }
        }
    }
}

impl NSKernel {
    /// Compute nonlinear RHS: P(ω × u) with dealiasing.
    /// Reads velocity from uhx/uhy/uhz, writes result to self.nl_x/nl_y/nl_z.
    fn nonlinear_rhs(
        &mut self,
        uhx: &[Complex64],
        uhy: &[Complex64],
        uhz: &[Complex64],
    ) {
        let n = self.n;
        let nzc = self.nzc;
        let slen = self.spec_len;
        let plen = self.phys_len;
        let norm = 1.0 / plen as f64;

        // 1. Compute vorticity in spectral space: ω̂ = ik × û
        for ix in 0..n {
            let kxi = self.kx[ix];
            for iy in 0..n {
                let kyi = self.ky[iy];
                for iz in 0..nzc {
                    let kzi = self.kz[iz];
                    let idx = ix * n * nzc + iy * nzc + iz;
                    let i = Complex64::new(0.0, 1.0);
                    self.omega_x[idx] = i * kyi * uhz[idx] - i * kzi * uhy[idx];
                    self.omega_y[idx] = i * kzi * uhx[idx] - i * kxi * uhz[idx];
                    self.omega_z[idx] = i * kxi * uhy[idx] - i * kyi * uhx[idx];
                }
            }
        }

        // 2. IFFT all 6 fields to physical space
        // phys_scratch[0..2] = velocity, phys_scratch[3..5] = vorticity
        unsafe {
            do_c2r(self.plan_c2r, uhx, &mut self.phys_scratch[0],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, uhy, &mut self.phys_scratch[1],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, uhz, &mut self.phys_scratch[2],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_x, &mut self.phys_scratch[3],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_y, &mut self.phys_scratch[4],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_z, &mut self.phys_scratch[5],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
        }

        // 3. Cross product in physical space: (ω × u)
        // Overwrite phys_scratch[3..5] with result
        for i in 0..plen {
            let ux = self.phys_scratch[0][i];
            let uy = self.phys_scratch[1][i];
            let uz = self.phys_scratch[2][i];
            let ox = self.phys_scratch[3][i];
            let oy = self.phys_scratch[4][i];
            let oz = self.phys_scratch[5][i];
            self.phys_scratch[3][i] = oy * uz - oz * uy;
            self.phys_scratch[4][i] = oz * ux - ox * uz;
            self.phys_scratch[5][i] = ox * uy - oy * ux;
        }

        // 4. FFT cross products to spectral space + dealias
        unsafe {
            do_r2c(self.plan_r2c, &self.phys_scratch[3], &mut self.nl_x,
                   &mut self.fft_phys, &mut self.fft_spec);
            do_r2c(self.plan_r2c, &self.phys_scratch[4], &mut self.nl_y,
                   &mut self.fft_phys, &mut self.fft_spec);
            do_r2c(self.plan_r2c, &self.phys_scratch[5], &mut self.nl_z,
                   &mut self.fft_phys, &mut self.fft_spec);
        }

        for i in 0..slen {
            let d = self.dealias[i];
            self.nl_x[i] *= d;
            self.nl_y[i] *= d;
            self.nl_z[i] *= d;
        }

        // 5. Leray projection: remove divergent part
        for ix in 0..n {
            let kxi = self.kx[ix];
            for iy in 0..n {
                let kyi = self.ky[iy];
                for iz in 0..nzc {
                    let kzi = self.kz[iz];
                    let idx = ix * n * nzc + iy * nzc + iz;
                    let k2s = self.k2_safe[idx];

                    let kdot = kxi * self.nl_x[idx]
                             + kyi * self.nl_y[idx]
                             + kzi * self.nl_z[idx];
                    let proj = kdot / k2s;

                    self.nl_x[idx] -= kxi * proj;
                    self.nl_y[idx] -= kyi * proj;
                    self.nl_z[idx] -= kzi * proj;
                }
            }
        }
    }
}

#[pymethods]
impl NSKernel {
    #[new]
    fn new(n: usize, nu: f64, dt: f64) -> PyResult<Self> {
        if n == 0 || n % 2 != 0 {
            return Err(pyo3::exceptions::PyValueError::new_err("n must be positive and even"));
        }

        let nzc = n / 2 + 1;
        let spec_len = n * n * nzc;
        let phys_len = n * n * n;

        // 1D wavenumber arrays
        let mut kx = vec![0.0; n];
        for i in 0..=n/2 { kx[i] = i as f64; }
        for i in (n/2+1)..n { kx[i] = (i as f64) - (n as f64); }
        let ky = kx.clone();
        let kz: Vec<f64> = (0..nzc).map(|i| i as f64).collect();

        // 3D precomputed arrays
        let mut k2 = vec![0.0; spec_len];
        let mut k2_safe = vec![0.0; spec_len];
        let mut dealias = vec![0.0; spec_len];
        let mut e_half = vec![0.0; spec_len];
        let mut e_full = vec![0.0; spec_len];
        let k_max = (n / 3) as f64;

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..nzc {
                    let idx = ix * n * nzc + iy * nzc + iz;
                    let ksq = kx[ix] * kx[ix] + ky[iy] * ky[iy] + kz[iz] * kz[iz];
                    k2[idx] = ksq;
                    k2_safe[idx] = if ksq == 0.0 { 1.0 } else { ksq };
                    dealias[idx] = if kx[ix].abs() <= k_max
                                   && ky[iy].abs() <= k_max
                                   && kz[iz].abs() <= k_max { 1.0 } else { 0.0 };
                    e_half[idx] = (-nu * ksq * dt * 0.5).exp();
                    e_full[idx] = (-nu * ksq * dt).exp();
                }
            }
        }

        // FFTW plans
        let mut plan_phys_buf = vec![0.0f64; phys_len];
        let mut plan_spec_buf = vec![Complex64::new(0.0, 0.0); spec_len];

        let plan_r2c;
        let plan_c2r;
        unsafe {
            plan_r2c = fftw_sys::fftw_plan_dft_r2c_3d(
                n as i32, n as i32, n as i32,
                plan_phys_buf.as_mut_ptr(),
                plan_spec_buf.as_mut_ptr(),
                fftw_sys::FFTW_PATIENT,
            );
            if plan_r2c.is_null() {
                return Err(pyo3::exceptions::PyRuntimeError::new_err("r2c plan creation failed"));
            }
            plan_c2r = fftw_sys::fftw_plan_dft_c2r_3d(
                n as i32, n as i32, n as i32,
                plan_spec_buf.as_mut_ptr(),
                plan_phys_buf.as_mut_ptr(),
                fftw_sys::FFTW_PATIENT,
            );
            if plan_c2r.is_null() {
                fftw_sys::fftw_destroy_plan(plan_r2c);
                return Err(pyo3::exceptions::PyRuntimeError::new_err("c2r plan creation failed"));
            }
        }

        let z = Complex64::new(0.0, 0.0);
        let phys_scratch = std::array::from_fn(|_| vec![0.0; phys_len]);

        Ok(NSKernel {
            n, nzc, spec_len, phys_len, nu, dt,
            t: 0.0, step_count: 0,
            kx, ky, kz, k2, k2_safe, dealias, e_half, e_full,
            ux_hat: vec![z; spec_len],
            uy_hat: vec![z; spec_len],
            uz_hat: vec![z; spec_len],
            plan_r2c, plan_c2r,
            _plan_phys_buf: plan_phys_buf,
            _plan_spec_buf: plan_spec_buf,
            fft_phys: vec![0.0; phys_len],
            fft_spec: vec![z; spec_len],
            phys_scratch,
            nl_x: vec![z; spec_len],
            nl_y: vec![z; spec_len],
            nl_z: vec![z; spec_len],
            omega_x: vec![z; spec_len],
            omega_y: vec![z; spec_len],
            omega_z: vec![z; spec_len],
            u0_x: vec![z; spec_len],
            u0_y: vec![z; spec_len],
            u0_z: vec![z; spec_len],
            acc_x: vec![z; spec_len],
            acc_y: vec![z; spec_len],
            acc_z: vec![z; spec_len],
        })
    }

    /// Set velocity from physical-space numpy arrays.
    fn set_velocity_physical<'py>(
        &mut self,
        ux: PyReadonlyArray3<'py, f64>,
        uy: PyReadonlyArray3<'py, f64>,
        uz: PyReadonlyArray3<'py, f64>,
    ) -> PyResult<()> {
        let shape = ux.shape();
        let n = self.n;
        if shape[0] != n || shape[1] != n || shape[2] != n {
            return Err(pyo3::exceptions::PyValueError::new_err(
                format!("expected shape ({n},{n},{n}), got ({},{},{})", shape[0], shape[1], shape[2])
            ));
        }

        let ux_s = ux.as_slice()?;
        let uy_s = uy.as_slice()?;
        let uz_s = uz.as_slice()?;

        unsafe {
            do_r2c(self.plan_r2c, ux_s, &mut self.ux_hat,
                   &mut self.fft_phys, &mut self.fft_spec);
            do_r2c(self.plan_r2c, uy_s, &mut self.uy_hat,
                   &mut self.fft_phys, &mut self.fft_spec);
            do_r2c(self.plan_r2c, uz_s, &mut self.uz_hat,
                   &mut self.fft_phys, &mut self.fft_spec);
        }

        self.t = 0.0;
        self.step_count = 0;
        Ok(())
    }

    /// One full IFRK4 step.
    fn step_rk4(&mut self) {
        let slen = self.spec_len;
        let dt = self.dt;

        // Save initial state
        self.u0_x.copy_from_slice(&self.ux_hat);
        self.u0_y.copy_from_slice(&self.uy_hat);
        self.u0_z.copy_from_slice(&self.uz_hat);

        // Stage 1: NL1 = NonlinearRHS(û₀)
        let uhx: Vec<Complex64> = self.ux_hat.clone();
        let uhy: Vec<Complex64> = self.uy_hat.clone();
        let uhz: Vec<Complex64> = self.uz_hat.clone();
        self.nonlinear_rhs(&uhx, &uhy, &uhz);

        // acc = (dt/6) * E_full * NL1
        // u1 = E_half * u0 + (dt/2) * E_half * NL1
        for i in 0..slen {
            let ef = self.e_full[i];
            let eh = self.e_half[i];
            self.acc_x[i] = (dt / 6.0 * ef) * self.nl_x[i];
            self.acc_y[i] = (dt / 6.0 * ef) * self.nl_y[i];
            self.acc_z[i] = (dt / 6.0 * ef) * self.nl_z[i];

            self.ux_hat[i] = eh * self.u0_x[i] + (dt * 0.5 * eh) * self.nl_x[i];
            self.uy_hat[i] = eh * self.u0_y[i] + (dt * 0.5 * eh) * self.nl_y[i];
            self.uz_hat[i] = eh * self.u0_z[i] + (dt * 0.5 * eh) * self.nl_z[i];
        }

        // Stage 2: NL2 = NonlinearRHS(u1)
        let uhx: Vec<Complex64> = self.ux_hat.clone();
        let uhy: Vec<Complex64> = self.uy_hat.clone();
        let uhz: Vec<Complex64> = self.uz_hat.clone();
        self.nonlinear_rhs(&uhx, &uhy, &uhz);

        // acc += (dt/3) * E_half * NL2
        // u2 = E_half * u0 + (dt/2) * E_half * NL2
        for i in 0..slen {
            let eh = self.e_half[i];
            self.acc_x[i] += (dt / 3.0 * eh) * self.nl_x[i];
            self.acc_y[i] += (dt / 3.0 * eh) * self.nl_y[i];
            self.acc_z[i] += (dt / 3.0 * eh) * self.nl_z[i];

            self.ux_hat[i] = eh * self.u0_x[i] + (dt * 0.5 * eh) * self.nl_x[i];
            self.uy_hat[i] = eh * self.u0_y[i] + (dt * 0.5 * eh) * self.nl_y[i];
            self.uz_hat[i] = eh * self.u0_z[i] + (dt * 0.5 * eh) * self.nl_z[i];
        }

        // Stage 3: NL3 = NonlinearRHS(u2)
        let uhx: Vec<Complex64> = self.ux_hat.clone();
        let uhy: Vec<Complex64> = self.uy_hat.clone();
        let uhz: Vec<Complex64> = self.uz_hat.clone();
        self.nonlinear_rhs(&uhx, &uhy, &uhz);

        // acc += (dt/3) * E_half * NL3
        // u3 = E_full * u0 + dt * E_half * NL3
        for i in 0..slen {
            let ef = self.e_full[i];
            let eh = self.e_half[i];
            self.acc_x[i] += (dt / 3.0 * eh) * self.nl_x[i];
            self.acc_y[i] += (dt / 3.0 * eh) * self.nl_y[i];
            self.acc_z[i] += (dt / 3.0 * eh) * self.nl_z[i];

            self.ux_hat[i] = ef * self.u0_x[i] + (dt * eh) * self.nl_x[i];
            self.uy_hat[i] = ef * self.u0_y[i] + (dt * eh) * self.nl_y[i];
            self.uz_hat[i] = ef * self.u0_z[i] + (dt * eh) * self.nl_z[i];
        }

        // Stage 4: NL4 = NonlinearRHS(u3)
        let uhx: Vec<Complex64> = self.ux_hat.clone();
        let uhy: Vec<Complex64> = self.uy_hat.clone();
        let uhz: Vec<Complex64> = self.uz_hat.clone();
        self.nonlinear_rhs(&uhx, &uhy, &uhz);

        // û_new = E_full * u0 + acc + (dt/6) * NL4
        for i in 0..slen {
            let ef = self.e_full[i];
            self.ux_hat[i] = ef * self.u0_x[i] + self.acc_x[i] + (dt / 6.0) * self.nl_x[i];
            self.uy_hat[i] = ef * self.u0_y[i] + self.acc_y[i] + (dt / 6.0) * self.nl_y[i];
            self.uz_hat[i] = ef * self.u0_z[i] + self.acc_z[i] + (dt / 6.0) * self.nl_z[i];
        }

        self.t += dt;
        self.step_count += 1;
    }

    /// Multiple RK4 steps.
    fn step_rk4_n(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step_rk4();
        }
    }

    /// Diagnostics via Parseval's theorem: (kinetic_energy, enstrophy).
    fn diagnostics(&self) -> (f64, f64) {
        let n = self.n;
        let nzc = self.nzc;
        let n3 = self.phys_len as f64;
        let mut energy = 0.0;
        let mut enstrophy = 0.0;

        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..nzc {
                    let idx = ix * n * nzc + iy * nzc + iz;
                    let u2 = self.ux_hat[idx].norm_sqr()
                           + self.uy_hat[idx].norm_sqr()
                           + self.uz_hat[idx].norm_sqr();

                    let factor = if iz == 0 || iz == n / 2 { 1.0 } else { 2.0 };
                    energy += factor * u2;
                    enstrophy += factor * self.k2[idx] * u2;
                }
            }
        }

        // Parseval (FFTW unnormalized):
        // sum_x |u(x)|² = (1/N³) sum_k |û(k)|²
        // E = (1/2) (1/N³) sum_x |u|² = (1/2) (1/N⁶) sum_k |û|²
        energy *= 0.5 / (n3 * n3);
        enstrophy *= 0.5 / (n3 * n3);
        (energy, enstrophy)
    }

    #[getter]
    fn t(&self) -> f64 { self.t }

    #[getter]
    fn step_count(&self) -> usize { self.step_count }

    /// Full diagnostics: E, Ω, Λ, K, |ω|_max, T1/T2 norms and ratios.
    /// Reuses scratch buffers; must only be called between timesteps.
    ///
    /// Returns a Python dict with keys:
    ///   E, Omega, Lambda, K, omega_max,
    ///   T1_inf, T1_L2, T2_inf, T2_L2, ratio_inf, ratio_L2
    fn diagnostics_full<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let n = self.n;
        let nzc = self.nzc;
        let slen = self.spec_len;
        let plen = self.phys_len;
        let norm = 1.0 / plen as f64;
        let eps = 1e-12;

        let (energy, enstrophy) = self.diagnostics();

        // --- Spectral vorticity: ω̂ = ik × û ---
        for ix in 0..n {
            let kxi = self.kx[ix];
            for iy in 0..n {
                let kyi = self.ky[iy];
                for iz in 0..nzc {
                    let kzi = self.kz[iz];
                    let idx = ix * n * nzc + iy * nzc + iz;
                    let ii = Complex64::new(0.0, 1.0);
                    self.omega_x[idx] = ii * kyi * self.uz_hat[idx] - ii * kzi * self.uy_hat[idx];
                    self.omega_y[idx] = ii * kzi * self.ux_hat[idx] - ii * kxi * self.uz_hat[idx];
                    self.omega_z[idx] = ii * kxi * self.uy_hat[idx] - ii * kyi * self.ux_hat[idx];
                }
            }
        }

        // --- IFFT vorticity → physical: phys_scratch[0..2] ---
        unsafe {
            do_c2r(self.plan_c2r, &self.omega_x, &mut self.phys_scratch[0],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_y, &mut self.phys_scratch[1],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_z, &mut self.phys_scratch[2],
                   &mut self.fft_phys, &mut self.fft_spec, norm);
        }

        // --- |ω| in phys_scratch[3], ω_max ---
        let mut omega_max = 0.0f64;
        for i in 0..plen {
            let ox = self.phys_scratch[0][i];
            let oy = self.phys_scratch[1][i];
            let oz = self.phys_scratch[2][i];
            let mag = (ox * ox + oy * oy + oz * oz).sqrt();
            self.phys_scratch[3][i] = mag;
            omega_max = omega_max.max(mag);
        }

        // --- ξ = ω / max(|ω|, ε): overwrite phys_scratch[0..2] ---
        for i in 0..plen {
            let mag = self.phys_scratch[3][i].max(eps);
            self.phys_scratch[0][i] /= mag;
            self.phys_scratch[1][i] /= mag;
            self.phys_scratch[2][i] /= mag;
        }

        // --- |∇ξ|²_F accumulator in phys_scratch[5] ---
        for i in 0..plen {
            self.phys_scratch[5][i] = 0.0;
        }

        for comp in 0..3usize {
            // FFT ξ_comp → omega_x (spectral scratch)
            unsafe {
                do_r2c(self.plan_r2c, &self.phys_scratch[comp], &mut self.omega_x,
                       &mut self.fft_phys, &mut self.fft_spec);
            }

            for dir in 0..3usize {
                // ik_dir × ξ̂_comp → nl_x
                for ix in 0..n {
                    for iy in 0..n {
                        for iz in 0..nzc {
                            let idx = ix * n * nzc + iy * nzc + iz;
                            let kj = match dir {
                                0 => self.kx[ix],
                                1 => self.ky[iy],
                                _ => self.kz[iz],
                            };
                            self.nl_x[idx] = Complex64::new(0.0, kj) * self.omega_x[idx];
                        }
                    }
                }

                // IFFT → phys_scratch[4]
                unsafe {
                    do_c2r(self.plan_c2r, &self.nl_x, &mut self.phys_scratch[4],
                           &mut self.fft_phys, &mut self.fft_spec, norm);
                }

                // Accumulate (∂ξ_comp/∂x_dir)²
                for i in 0..plen {
                    let v = self.phys_scratch[4][i];
                    self.phys_scratch[5][i] += v * v;
                }
            }
        }

        // --- Λ = sup |∇ξ|_F ---
        let mut lambda = 0.0f64;
        for i in 0..plen {
            lambda = lambda.max(self.phys_scratch[5][i].sqrt());
        }

        // --- η = (ν³/ε_diss)^(1/4), K = Λ × η ---
        let eps_diss = 2.0 * self.nu * enstrophy;
        let eta = if eps_diss > eps {
            (self.nu.powi(3) / eps_diss).powf(0.25)
        } else {
            f64::INFINITY
        };
        let k_metric = lambda * eta;

        // --- ∇|ω| for T1: FFT |ω| → omega_x, gradient accumulator in phys_scratch[0] ---
        // (ξ no longer needed in phys_scratch[0..2])
        unsafe {
            do_r2c(self.plan_r2c, &self.phys_scratch[3], &mut self.omega_x,
                   &mut self.fft_phys, &mut self.fft_spec);
        }

        for i in 0..plen {
            self.phys_scratch[0][i] = 0.0; // |∇|ω||² accumulator
        }

        for dir in 0..3usize {
            for ix in 0..n {
                for iy in 0..n {
                    for iz in 0..nzc {
                        let idx = ix * n * nzc + iy * nzc + iz;
                        let kj = match dir {
                            0 => self.kx[ix],
                            1 => self.ky[iy],
                            _ => self.kz[iz],
                        };
                        self.nl_x[idx] = Complex64::new(0.0, kj) * self.omega_x[idx];
                    }
                }
            }

            unsafe {
                do_c2r(self.plan_c2r, &self.nl_x, &mut self.phys_scratch[4],
                       &mut self.fft_phys, &mut self.fft_spec, norm);
            }

            for i in 0..plen {
                let v = self.phys_scratch[4][i];
                self.phys_scratch[0][i] += v * v;
            }
        }

        // --- T1/T2 norms ---
        // |T1|_F = |∇|ω|| (since |ξ|=1, outer product norm = |a|·|b|)
        // |T2|_F = |ω| × |∇ξ|_F
        let mut t1_inf = 0.0f64;
        let mut t1_l2_sq = 0.0;
        let mut t2_inf = 0.0f64;
        let mut t2_l2_sq = 0.0;

        for i in 0..plen {
            let t1 = self.phys_scratch[0][i].sqrt();
            let t2 = self.phys_scratch[3][i] * self.phys_scratch[5][i].sqrt();
            t1_inf = t1_inf.max(t1);
            t1_l2_sq += t1 * t1;
            t2_inf = t2_inf.max(t2);
            t2_l2_sq += t2 * t2;
        }

        let t1_l2 = (t1_l2_sq / plen as f64).sqrt();
        let t2_l2 = (t2_l2_sq / plen as f64).sqrt();
        let ratio_inf = t2_inf / t1_inf.max(eps);
        let ratio_l2 = t2_l2 / t1_l2.max(eps);

        // --- Pack into Python dict ---
        let dict = PyDict::new(py);
        dict.set_item("E", energy)?;
        dict.set_item("Omega", enstrophy)?;
        dict.set_item("Lambda", lambda)?;
        dict.set_item("K", k_metric)?;
        dict.set_item("omega_max", omega_max)?;
        dict.set_item("T1_inf", t1_inf)?;
        dict.set_item("T1_L2", t1_l2)?;
        dict.set_item("T2_inf", t2_inf)?;
        dict.set_item("T2_L2", t2_l2)?;
        dict.set_item("ratio_inf", ratio_inf)?;
        dict.set_item("ratio_L2", ratio_l2)?;
        Ok(dict)
    }

    /// Get snapshot fields as numpy arrays for HDF5 export.
    /// Returns dict with velocity(3), vorticity(3), omega_mag(1), xi(3),
    /// T1_norm_field(1), T2_norm_field(1) — all shape (N,N,N).
    fn get_snapshot_fields<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let n = self.n;
        let nzc = self.nzc;
        let slen = self.spec_len;
        let plen = self.phys_len;
        let norm = 1.0 / plen as f64;
        let eps = 1e-12;

        // --- Velocity in physical space ---
        let mut vel_x = vec![0.0f64; plen];
        let mut vel_y = vec![0.0f64; plen];
        let mut vel_z = vec![0.0f64; plen];
        unsafe {
            do_c2r(self.plan_c2r, &self.ux_hat, &mut vel_x,
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.uy_hat, &mut vel_y,
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.uz_hat, &mut vel_z,
                   &mut self.fft_phys, &mut self.fft_spec, norm);
        }

        // --- Spectral vorticity → physical ---
        for ix in 0..n {
            let kxi = self.kx[ix];
            for iy in 0..n {
                let kyi = self.ky[iy];
                for iz in 0..nzc {
                    let kzi = self.kz[iz];
                    let idx = ix * n * nzc + iy * nzc + iz;
                    let ii = Complex64::new(0.0, 1.0);
                    self.omega_x[idx] = ii * kyi * self.uz_hat[idx] - ii * kzi * self.uy_hat[idx];
                    self.omega_y[idx] = ii * kzi * self.ux_hat[idx] - ii * kxi * self.uz_hat[idx];
                    self.omega_z[idx] = ii * kxi * self.uy_hat[idx] - ii * kyi * self.ux_hat[idx];
                }
            }
        }

        let mut vort_x = vec![0.0f64; plen];
        let mut vort_y = vec![0.0f64; plen];
        let mut vort_z = vec![0.0f64; plen];
        unsafe {
            do_c2r(self.plan_c2r, &self.omega_x, &mut vort_x,
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_y, &mut vort_y,
                   &mut self.fft_phys, &mut self.fft_spec, norm);
            do_c2r(self.plan_c2r, &self.omega_z, &mut vort_z,
                   &mut self.fft_phys, &mut self.fft_spec, norm);
        }

        // --- |ω|, ξ ---
        let mut omega_mag = vec![0.0f64; plen];
        let mut xi_x = vec![0.0f64; plen];
        let mut xi_y = vec![0.0f64; plen];
        let mut xi_z = vec![0.0f64; plen];
        for i in 0..plen {
            let mag = (vort_x[i] * vort_x[i] + vort_y[i] * vort_y[i] + vort_z[i] * vort_z[i]).sqrt();
            omega_mag[i] = mag;
            let denom = mag.max(eps);
            xi_x[i] = vort_x[i] / denom;
            xi_y[i] = vort_y[i] / denom;
            xi_z[i] = vort_z[i] / denom;
        }

        // --- |∇ξ|²_F ---
        let mut grad_xi_sq = vec![0.0f64; plen];
        for comp_data in [&xi_x, &xi_y, &xi_z] {
            unsafe {
                do_r2c(self.plan_r2c, comp_data, &mut self.omega_x,
                       &mut self.fft_phys, &mut self.fft_spec);
            }
            for dir in 0..3usize {
                for ix in 0..n {
                    for iy in 0..n {
                        for iz in 0..nzc {
                            let idx = ix * n * nzc + iy * nzc + iz;
                            let kj = match dir {
                                0 => self.kx[ix],
                                1 => self.ky[iy],
                                _ => self.kz[iz],
                            };
                            self.nl_x[idx] = Complex64::new(0.0, kj) * self.omega_x[idx];
                        }
                    }
                }
                unsafe {
                    do_c2r(self.plan_c2r, &self.nl_x, &mut self.phys_scratch[4],
                           &mut self.fft_phys, &mut self.fft_spec, norm);
                }
                for i in 0..plen {
                    let v = self.phys_scratch[4][i];
                    grad_xi_sq[i] += v * v;
                }
            }
        }

        // --- ∇|ω| for T1 ---
        unsafe {
            do_r2c(self.plan_r2c, &omega_mag, &mut self.omega_x,
                   &mut self.fft_phys, &mut self.fft_spec);
        }
        let mut grad_omag_sq = vec![0.0f64; plen];
        for dir in 0..3usize {
            for ix in 0..n {
                for iy in 0..n {
                    for iz in 0..nzc {
                        let idx = ix * n * nzc + iy * nzc + iz;
                        let kj = match dir {
                            0 => self.kx[ix],
                            1 => self.ky[iy],
                            _ => self.kz[iz],
                        };
                        self.nl_x[idx] = Complex64::new(0.0, kj) * self.omega_x[idx];
                    }
                }
            }
            unsafe {
                do_c2r(self.plan_c2r, &self.nl_x, &mut self.phys_scratch[4],
                       &mut self.fft_phys, &mut self.fft_spec, norm);
            }
            for i in 0..plen {
                let v = self.phys_scratch[4][i];
                grad_omag_sq[i] += v * v;
            }
        }

        // --- T1_norm = |∇|ω||, T2_norm = |ω| × |∇ξ|_F ---
        let mut t1_norm = vec![0.0f64; plen];
        let mut t2_norm = vec![0.0f64; plen];
        for i in 0..plen {
            t1_norm[i] = grad_omag_sq[i].sqrt();
            t2_norm[i] = omega_mag[i] * grad_xi_sq[i].sqrt();
        }

        // --- Pack into dict of numpy arrays ---
        let dict = PyDict::new(py);
        let shape = [n, n, n];

        macro_rules! add_field {
            ($name:expr, $data:expr) => {
                let arr = $data.into_pyarray(py).reshape(shape)
                    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
                dict.set_item($name, arr)?;
            };
        }

        add_field!("ux", vel_x);
        add_field!("uy", vel_y);
        add_field!("uz", vel_z);
        add_field!("omega_x", vort_x);
        add_field!("omega_y", vort_y);
        add_field!("omega_z", vort_z);
        add_field!("omega_mag", omega_mag);
        add_field!("xi_x", xi_x);
        add_field!("xi_y", xi_y);
        add_field!("xi_z", xi_z);
        add_field!("T1_norm", t1_norm);
        add_field!("T2_norm", t2_norm);

        Ok(dict)
    }
}
