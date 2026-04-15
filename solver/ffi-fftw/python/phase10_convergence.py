#!/usr/bin/env python3
"""Phase 10b: Convergence runs at N=64, N=128, N=256.

Unattended overnight script with thermal/disk/wall-clock monitoring.
Antiparallel vortex tubes IC, Re=6400, full Λ/K/T1/T2 diagnostics.

Usage:
    python phase10_convergence.py
"""
import sys
import os
import csv
import json
import time
import glob
import resource
import subprocess

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import h5py
import fftw_training
from fftw_training import NSKernel
from antiparallel_tubes_ic import antiparallel_tubes_ic

# --- Constants ---
NU = 2 * np.pi / 6400
DT = 0.001
DIAG_DT = 0.01  # diagnostics every 0.01 sim time
STEPS_PER_DIAG = round(DIAG_DT / DT)  # = 10

# Per-resolution config (handoff + patches)
RUN_CONFIG = {
    64:  {"t_max": 40.0, "snapshots": [0.0, 10.0, 20.0, 30.0, 36.0, 40.0]},
    128: {"t_max": 40.0, "snapshots": [0.0, 10.0, 20.0, 30.0, 36.0, 40.0]},
    256: {"t_max": 30.0, "snapshots": [0.0, 10.0, 20.0, 25.0, 28.0, 30.0]},
}

# Probe wall-clock per step (from probe runs)
PROBE_WALL_PER_STEP = {64: 0.0454, 128: 0.3885, 256: 3.8784}

# Monitoring thresholds (patched)
DISK_PRE_RUN_GB = 15.0
DISK_POST_SNAP_GB = 8.0
THERMAL_LIMIT_C = 90.0
THERMAL_SUSTAIN_S = 60.0
WALL_CLOCK_FACTOR = 1.5
WALL_CLOCK_CONSEC = 20

SNAPSHOT_DIR = "snapshots"
DIAG_FIELDS = [
    "t", "E", "Omega", "Lambda", "K", "omega_max",
    "T1_inf", "T1_L2", "T2_inf", "T2_L2", "ratio_inf", "ratio_L2",
]

# Reporting times for LOG.md table
REPORT_TIMES = [5, 10, 14, 20, 25, 30, 35, 36, 38, 40]


def get_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def get_disk_free_gb():
    st = os.statvfs(".")
    return (st.f_bavail * st.f_frsize) / (1024**3)


def get_cpu_temp_c():
    """Read CPU package temp from sysfs. Returns None if unavailable."""
    for zone_dir in sorted(glob.glob("/sys/class/thermal/thermal_zone*/")):
        try:
            with open(os.path.join(zone_dir, "type")) as f:
                ztype = f.read().strip()
            if ztype == "x86_pkg_temp":
                with open(os.path.join(zone_dir, "temp")) as f:
                    return int(f.read().strip()) / 1000.0
        except (IOError, ValueError):
            continue
    return None


def get_swap_activity():
    pswpin = pswpout = 0
    with open("/proc/vmstat") as f:
        for line in f:
            if line.startswith("pswpin"):
                pswpin = int(line.split()[1])
            elif line.startswith("pswpout"):
                pswpout = int(line.split()[1])
    return pswpin, pswpout


def write_failure(n, step, t, last_diag, error_msg, wall_clock, peak_rss):
    path = f"failure_{n}.json"
    data = {
        "N": n,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "last_step": step,
        "last_t": t,
        "last_diagnostics": last_diag,
        "error": error_msg,
        "peak_rss_gb": round(peak_rss, 3),
        "wall_clock_s": round(wall_clock, 1),
        "disk_free_gb": round(get_disk_free_gb(), 1),
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  Wrote {path}")


def run_single(n):
    """Run convergence for a single resolution. Returns True on success."""
    cfg = RUN_CONFIG[n]
    t_max = cfg["t_max"]
    snap_times = set(cfg["snapshots"])
    total_steps = int(t_max / DT)
    total_diags = int(t_max / DIAG_DT) + 1
    probe_wps = PROBE_WALL_PER_STEP.get(n, 999)

    print(f"\n{'='*70}")
    print(f"CONVERGENCE RUN N={n}, t_max={t_max}, dt={DT}, {total_steps} steps")
    print(f"{'='*70}")

    # --- Disk gate ---
    disk_free = get_disk_free_gb()
    print(f"Disk free: {disk_free:.1f} GB (need ≥{DISK_PRE_RUN_GB})")
    if disk_free < DISK_PRE_RUN_GB:
        print(f"ABORT: insufficient disk ({disk_free:.1f} < {DISK_PRE_RUN_GB})")
        return False

    # --- Create kernel ---
    fftw_training.load_wisdom("wisdom.dat")
    t0_plan = time.monotonic()
    kernel = NSKernel(n, NU, DT)
    plan_time = time.monotonic() - t0_plan
    fftw_training.save_wisdom("wisdom.dat")
    print(f"Kernel created in {plan_time:.1f}s")

    # --- Set IC ---
    ux, uy, uz = antiparallel_tubes_ic(n)
    kernel.set_velocity_physical(ux, uy, uz)
    del ux, uy, uz

    peak_rss = get_rss_gb()
    disk_start = get_disk_free_gb()
    print(f"Baseline RSS: {peak_rss:.3f} GB")

    # --- Open CSV ---
    csv_path = f"run_{n}.csv"
    csv_file = open(csv_path, "w", newline="")
    csv_writer = csv.DictWriter(csv_file, fieldnames=DIAG_FIELDS)
    csv_writer.writeheader()

    # --- Tracking ---
    os.makedirs(SNAPSHOT_DIR, exist_ok=True)
    hdf5_bytes_total = 0
    temps = []
    last_diag = None
    slow_step_count = 0
    thermal_over_start = None
    snap_done = set()

    t_run_start = time.monotonic()

    def maybe_snapshot(sim_t):
        """Write HDF5 snapshot if sim_t matches a snapshot time."""
        nonlocal hdf5_bytes_total
        for st in snap_times:
            if st in snap_done:
                continue
            if abs(sim_t - st) < DT * 0.5:
                snap_path = f"{SNAPSHOT_DIR}/run_{n}_t{st:06.1f}.h5"
                print(f"  Snapshot at t={st:.1f} → {snap_path}")
                fields = kernel.get_snapshot_fields()
                with h5py.File(snap_path, "w") as hf:
                    hf.attrs["N"] = n
                    hf.attrs["t"] = sim_t
                    hf.attrs["nu"] = NU
                    hf.attrs["dt"] = DT
                    for key, arr in fields.items():
                        hf.create_dataset(key, data=np.array(arr), compression="gzip",
                                          compression_opts=4)
                snap_size = os.path.getsize(snap_path)
                hdf5_bytes_total += snap_size
                print(f"    {snap_size / 1024**2:.1f} MB")
                snap_done.add(st)

                # Post-snapshot disk check
                df = get_disk_free_gb()
                if df < DISK_POST_SNAP_GB:
                    raise RuntimeError(f"Disk free {df:.1f} GB < {DISK_POST_SNAP_GB} after snapshot")
                return

    try:
        for diag_idx in range(total_diags):
            sim_t = kernel.t

            # --- Full diagnostics ---
            d = kernel.diagnostics_full()
            row = {"t": sim_t}
            for k in DIAG_FIELDS[1:]:
                row[k] = d[k]
            csv_writer.writerow(row)
            last_diag = dict(row)

            # --- Snapshot check ---
            maybe_snapshot(sim_t)

            # --- Progress ---
            if diag_idx % 100 == 0:
                elapsed = time.monotonic() - t_run_start
                temp = get_cpu_temp_c()
                rss = get_rss_gb()
                peak_rss = max(peak_rss, rss)
                temp_str = f"{temp:.0f}°C" if temp else "N/A"
                print(f"  t={sim_t:7.2f}  E={d['E']:.8f}  Λ={d['Lambda']:10.2f}  "
                      f"K={d['K']:.4f}  RSS={rss:.2f}GB  T={temp_str}  "
                      f"wall={elapsed:.0f}s")
                if temp is not None:
                    temps.append(temp)

            # --- Halt checks (every diagnostic interval) ---
            # Thermal
            temp = get_cpu_temp_c()
            if temp is not None:
                temps.append(temp)
                if temp > THERMAL_LIMIT_C:
                    if thermal_over_start is None:
                        thermal_over_start = time.monotonic()
                    elif time.monotonic() - thermal_over_start > THERMAL_SUSTAIN_S:
                        raise RuntimeError(
                            f"CPU temp {temp:.0f}°C > {THERMAL_LIMIT_C}°C "
                            f"sustained for {THERMAL_SUSTAIN_S}s")
                else:
                    thermal_over_start = None

            # NaN check
            if not np.isfinite(d["E"]) or not np.isfinite(d["Omega"]):
                raise RuntimeError(f"NaN/Inf in diagnostics at t={sim_t}")
            if d["E"] < 0:
                raise RuntimeError(f"Negative energy at t={sim_t}: E={d['E']}")

            # Done?
            if sim_t >= t_max - DT * 0.5:
                break

            # --- Step ---
            t_step_start = time.monotonic()
            kernel.step_rk4_n(STEPS_PER_DIAG)
            step_wall = (time.monotonic() - t_step_start) / STEPS_PER_DIAG

            # Wall-clock halt
            if step_wall > WALL_CLOCK_FACTOR * probe_wps:
                slow_step_count += 1
                if slow_step_count >= WALL_CLOCK_CONSEC:
                    raise RuntimeError(
                        f"Wall-clock {step_wall:.4f}s/step > "
                        f"{WALL_CLOCK_FACTOR}× probe ({probe_wps:.4f}s) "
                        f"for {WALL_CLOCK_CONSEC} consecutive batches")
            else:
                slow_step_count = 0

        # --- Final snapshot if not yet done ---
        maybe_snapshot(kernel.t)

    except Exception as e:
        wall_total = time.monotonic() - t_run_start
        print(f"\n  HALT: {e}")
        write_failure(n, kernel.step_count, kernel.t, last_diag, str(e),
                      wall_total, peak_rss)
        csv_file.close()

        # Check if this is a kernel bug (panic, compile error) vs isolated failure
        if "panic" in str(e).lower() or "symbol" in str(e).lower():
            print("  KERNEL BUG DETECTED — halting all remaining resolutions")
            return "kernel_bug"
        return False

    wall_total = time.monotonic() - t_run_start
    csv_file.close()

    disk_end = get_disk_free_gb()
    temp_min = min(temps) if temps else 0
    temp_max = max(temps) if temps else 0
    temp_mean = np.mean(temps) if temps else 0

    print(f"\n  COMPLETED: N={n}, t={kernel.t:.2f}, {kernel.step_count} steps")
    print(f"  Wall clock: {wall_total:.1f}s ({wall_total/3600:.2f}h)")
    print(f"  Peak RSS: {peak_rss:.3f} GB")
    print(f"  CPU temp: min={temp_min:.0f} max={temp_max:.0f} mean={temp_mean:.0f}°C")
    print(f"  Disk: start={disk_start:.1f} end={disk_end:.1f} GB")
    print(f"  HDF5 total: {hdf5_bytes_total / 1024**2:.1f} MB")
    print(f"  Wrote {csv_path}")

    # --- Summary JSON ---
    summary = {
        "N": n, "t_max": t_max, "t_final": kernel.t,
        "steps": kernel.step_count, "wall_clock_s": round(wall_total, 1),
        "peak_rss_gb": round(peak_rss, 3),
        "temp_min_c": round(temp_min, 1), "temp_max_c": round(temp_max, 1),
        "temp_mean_c": round(temp_mean, 1),
        "disk_start_gb": round(disk_start, 1), "disk_end_gb": round(disk_end, 1),
        "hdf5_bytes": hdf5_bytes_total,
    }
    with open(f"run_{n}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    return True


def plot_diagnostics(n, csv_path):
    """Generate diagnostic plots from CSV."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t, E, Om, Lam, K, om_max = [], [], [], [], [], []
    ratio_inf, ratio_l2 = [], []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            t.append(float(row["t"]))
            E.append(float(row["E"]))
            Om.append(float(row["Omega"]))
            Lam.append(float(row["Lambda"]))
            K.append(float(row["K"]))
            om_max.append(float(row["omega_max"]))
            ratio_inf.append(float(row["ratio_inf"]))
            ratio_l2.append(float(row["ratio_L2"]))

    t = np.array(t)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0, 0].plot(t, Lam, "b-", lw=1)
    axes[0, 0].set_ylabel("Λ = sup|∇ξ|_F")
    axes[0, 0].set_title(f"Λ(t) — N={n}")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(t, K, "r-", lw=1)
    axes[0, 1].set_ylabel("K = Λ × η")
    axes[0, 1].set_title(f"K(t) — N={n}")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(t, Om, "k-", lw=1)
    axes[1, 0].set_ylabel("Enstrophy Ω")
    axes[1, 0].set_title(f"Enstrophy(t) — N={n}")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(t, ratio_inf, "b-", lw=1, label="||T2||_∞/||T1||_∞")
    axes[1, 1].plot(t, ratio_l2, "r--", lw=1, label="||T2||_L²/||T1||_L²")
    axes[1, 1].set_ylabel("Ratio")
    axes[1, 1].set_title(f"||T2||/||T1|| — N={n}")
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    for ax in axes.flat:
        ax.set_xlabel("t")

    plt.tight_layout()
    out = f"run_{n}_diagnostics.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Wrote {out}")


def extract_report_values(csv_path, t_max):
    """Extract diagnostic values at reporting times for LOG.md table."""
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})

    report = {}
    for rt in REPORT_TIMES:
        if rt > t_max + 0.1:
            continue
        # Find closest row
        best_idx = min(range(len(rows)), key=lambda i: abs(rows[i]["t"] - rt))
        if abs(rows[best_idx]["t"] - rt) < 0.1:
            report[rt] = rows[best_idx]
    return report


if __name__ == "__main__":
    print(f"Phase 10b Convergence Runs — {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"ν = {NU:.6e}, dt = {DT}")
    print(f"Disk free: {get_disk_free_gb():.1f} GB")
    temp = get_cpu_temp_c()
    print(f"CPU temp: {temp:.0f}°C" if temp else "CPU temp: N/A")
    print()

    results = {}
    run_order = [64, 128, 256]

    for n in run_order:
        result = run_single(n)

        if result == "kernel_bug":
            print(f"\nKERNEL BUG at N={n} — halting all remaining resolutions.")
            break

        if result:
            plot_diagnostics(n, f"run_{n}.csv")
            results[n] = "completed"
            # Commit after each successful run
            os.system(f'git add run_{n}.csv run_{n}_summary.json run_{n}_diagnostics.png '
                      f'snapshots/run_{n}_*.h5 2>/dev/null; '
                      f'git commit -m "phase10b:convergence_N{n}\n\n'
                      f'Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>"')
        else:
            results[n] = "failed"
            # Still commit failure forensics
            os.system(f'git add failure_{n}.json run_{n}.csv 2>/dev/null; '
                      f'git commit -m "phase10b:convergence_N{n}_failed\n\n'
                      f'Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>"')
            # Continue to next resolution (isolated failure)

        # Free kernel memory before next resolution
        import gc
        gc.collect()

    # --- Final status report ---
    print(f"\n{'='*70}")
    print("PHASE 10b STATUS REPORT")
    print(f"{'='*70}")
    for n in run_order:
        status = results.get(n, "not_started")
        print(f"  N={n}: {status}")
    print(f"\nTotal disk usage: {get_disk_free_gb():.1f} GB free remaining")
    print(f"Git log:")
    os.system("git log --oneline -10")

    if all(v == "completed" for v in results.values()):
        print("\nPhase 10b complete, awaiting offline analysis handoff.")
    else:
        failed = [n for n, v in results.items() if v != "completed"]
        print(f"\nPhase 10b halted/failed at: {failed}")
        print(f"Forensics in: {['failure_'+str(n)+'.json' for n in failed]}")
