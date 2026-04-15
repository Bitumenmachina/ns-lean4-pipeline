mod kernel;

use numpy::{PyArray1, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray3, PyUntypedArrayMethods, IntoPyArray};
use num_complex::Complex64;
use pyo3::prelude::*;
use std::ffi::CString;

#[pyfunction]
fn rfft_1d<'py>(py: Python<'py>, input: PyReadonlyArray1<'py, f64>) -> PyResult<Bound<'py, PyArray1<Complex64>>> {
    let slice = input.as_slice()?;
    let n = slice.len();
    let n_complex = n / 2 + 1;

    let mut in_buf: Vec<f64> = slice.to_vec();
    let mut out_buf: Vec<fftw_sys::fftw_complex> = vec![Complex64::new(0.0, 0.0); n_complex];

    unsafe {
        let plan = fftw_sys::fftw_plan_dft_r2c_1d(
            n as i32,
            in_buf.as_mut_ptr(),
            out_buf.as_mut_ptr(),
            fftw_sys::FFTW_ESTIMATE,
        );
        if plan.is_null() {
            return Err(pyo3::exceptions::PyRuntimeError::new_err("FFTW plan creation failed"));
        }
        fftw_sys::fftw_execute(plan);
        fftw_sys::fftw_destroy_plan(plan);
    }

    Ok(out_buf.into_pyarray(py))
}

#[pyfunction]
#[pyo3(signature = (input, patient=true))]
fn rfft_3d<'py>(
    py: Python<'py>,
    input: PyReadonlyArray3<'py, f64>,
    patient: bool,
) -> PyResult<Bound<'py, PyArray3<Complex64>>> {
    let shape = input.shape();
    let nx = shape[0];
    let ny = shape[1];
    let nz = shape[2];
    let nz_complex = nz / 2 + 1;

    let in_len = nx * ny * nz;
    let out_len = nx * ny * nz_complex;

    // Allocate buffers for planning (FFTW_PATIENT overwrites buffers during planning)
    let mut in_buf: Vec<f64> = vec![0.0; in_len];
    let mut out_buf: Vec<fftw_sys::fftw_complex> = vec![Complex64::new(0.0, 0.0); out_len];

    let flags = if patient {
        fftw_sys::FFTW_PATIENT
    } else {
        fftw_sys::FFTW_ESTIMATE
    };

    unsafe {
        let plan = fftw_sys::fftw_plan_dft_r2c_3d(
            nx as i32,
            ny as i32,
            nz as i32,
            in_buf.as_mut_ptr(),
            out_buf.as_mut_ptr(),
            flags,
        );
        if plan.is_null() {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "FFTW 3D plan creation failed",
            ));
        }

        // Copy input data AFTER planning (PATIENT/MEASURE overwrite buffers)
        let slice = input.as_slice()?;
        in_buf.copy_from_slice(slice);

        fftw_sys::fftw_execute(plan);
        fftw_sys::fftw_destroy_plan(plan);
    }

    let flat = out_buf.into_pyarray(py);
    let reshaped = flat.reshape([nx, ny, nz_complex])
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("reshape failed: {e}")))?;
    Ok(reshaped.into())
}

#[pyfunction]
#[pyo3(signature = (input, nz_full))]
fn irfft_3d<'py>(
    py: Python<'py>,
    input: PyReadonlyArray3<'py, Complex64>,
    nz_full: usize,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    let shape = input.shape();
    let nx = shape[0];
    let ny = shape[1];
    let nz_complex = shape[2];

    if nz_complex != nz_full / 2 + 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "input shape[2]={} inconsistent with nz_full={} (expected {})",
            nz_complex, nz_full, nz_full / 2 + 1
        )));
    }

    let in_len = nx * ny * nz_complex;
    let out_len = nx * ny * nz_full;

    // Allocate buffers for planning (c2r destroys input during both planning and execution)
    let mut in_buf: Vec<fftw_sys::fftw_complex> = vec![Complex64::new(0.0, 0.0); in_len];
    let mut out_buf: Vec<f64> = vec![0.0; out_len];

    unsafe {
        let plan = fftw_sys::fftw_plan_dft_c2r_3d(
            nx as i32,
            ny as i32,
            nz_full as i32,
            in_buf.as_mut_ptr(),
            out_buf.as_mut_ptr(),
            fftw_sys::FFTW_ESTIMATE,
        );
        if plan.is_null() {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "FFTW 3D c2r plan creation failed",
            ));
        }

        // Copy input data AFTER planning (c2r destroys input buffer)
        let slice = input.as_slice()?;
        in_buf.copy_from_slice(slice);

        fftw_sys::fftw_execute(plan);
        fftw_sys::fftw_destroy_plan(plan);
    }

    // FFTW c2r is unnormalized: output = N³ * true_inverse
    let n3 = (nx * ny * nz_full) as f64;
    for v in out_buf.iter_mut() {
        *v /= n3;
    }

    let flat = out_buf.into_pyarray(py);
    let reshaped = flat.reshape([nx, ny, nz_full])
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("reshape failed: {e}")))?;
    Ok(reshaped.into())
}

#[pyfunction]
fn save_wisdom(path: &str) -> PyResult<bool> {
    let c_path = CString::new(path)
        .map_err(|_| pyo3::exceptions::PyValueError::new_err("invalid path"))?;
    unsafe {
        let ret = fftw_sys::fftw_export_wisdom_to_filename(c_path.as_ptr());
        Ok(ret != 0)
    }
}

#[pyfunction]
fn load_wisdom(path: &str) -> PyResult<bool> {
    let c_path = CString::new(path)
        .map_err(|_| pyo3::exceptions::PyValueError::new_err("invalid path"))?;
    unsafe {
        let ret = fftw_sys::fftw_import_wisdom_from_filename(c_path.as_ptr());
        Ok(ret != 0)
    }
}

#[pymodule]
fn fftw_training(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rfft_1d, m)?)?;
    m.add_function(wrap_pyfunction!(rfft_3d, m)?)?;
    m.add_function(wrap_pyfunction!(irfft_3d, m)?)?;
    m.add_function(wrap_pyfunction!(save_wisdom, m)?)?;
    m.add_function(wrap_pyfunction!(load_wisdom, m)?)?;
    m.add_class::<kernel::NSKernel>()?;
    Ok(())
}
