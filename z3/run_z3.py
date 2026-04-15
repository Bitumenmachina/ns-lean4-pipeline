#!/usr/bin/env python3
"""Run Z3 on all certificates and collect results."""
import json, subprocess, sys
from pathlib import Path

Z3_DIR = Path(__file__).parent
smt2 = Z3_DIR / "all_certs.smt2"

try:
    import z3
    HAS_Z3_PY = True
except ImportError:
    HAS_Z3_PY = False

cert_names = ['K_t40_gt_1', 'K_crossover', 'convergence_ratio_bounded', 'convergence_ratios_decreasing', 'onset_ft_strictly_increasing', 'staircase_peaks_increasing', 'Lambda_grows_Omega_decays', 'subcritical_tame', 'K_def_consistent', 'regime_boundaries_ordered']

def run_z3_cli():
    r = subprocess.run(["z3", str(smt2)], capture_output=True, text=True)
    lines = [l.strip() for l in r.stdout.strip().split("\n") if l.strip()]
    return lines

def run_z3_python():
    """Parse SMT2 and use z3 Python bindings."""
    import z3
    text = smt2.read_text()
    results = []

    # Split by push/pop blocks
    blocks = text.split("(push)")
    for block in blocks[1:]:  # skip header
        inner = block.split("(pop)")[0]
        s = z3.Solver()
        s.set("timeout", 10000)
        s.from_string("(set-logic QF_LRA)\n" + inner)
        r = s.check()
        results.append("sat" if r == z3.sat else "unsat" if r == z3.unsat else "unknown")
    return results

if HAS_Z3_PY:
    lines = run_z3_python()
else:
    lines = run_z3_cli()

results = {}
for i, name in enumerate(cert_names):
    status = lines[i] if i < len(lines) else "error"
    results[name] = {"status": status, "index": i}
    mark = "✓" if status == "sat" else "✗"
    print(f"  {mark} {name}: {status}")

out = Z3_DIR / "results.json"
out.write_text(json.dumps(results, indent=2))
print(f"\nWrote {out}")

failed = [n for n, r in results.items() if r["status"] != "sat"]
if failed:
    print(f"\nFAILED: {failed}")
    sys.exit(1)
else:
    print(f"\nAll {len(results)} certificates SAT")
