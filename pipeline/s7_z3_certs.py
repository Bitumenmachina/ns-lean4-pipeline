#!/usr/bin/env python3
"""S7: Generate Z3 SMT-LIB2 certificate scripts from params.json.

Each assertion encodes a numerical fact from the V5 experiments.
Z3 proves these in QF_LRA (quantifier-free linear real arithmetic).
Results feed into S8 to become Lean axioms.

Output: z3/all_certs.smt2, z3/results.json
"""

import json
from pathlib import Path

PARAMS = Path(__file__).parent.parent / "params.json"
Z3_DIR = Path(__file__).parent.parent / "z3"


def main():
    p = json.load(PARAMS.open())
    Z3_DIR.mkdir(parents=True, exist_ok=True)

    blocks = []

    # ── Header ───────────────────────────────────────────────────────
    blocks.append("; NS-Lean4 Z3 Certificates")
    blocks.append("; Auto-generated from params.json")
    blocks.append("; Logic: QF_LRA (quantifier-free linear real arithmetic)")
    blocks.append("(set-logic QF_LRA)")
    blocks.append("")

    certs = []  # (name, description, smt2_block)

    # ── CERT-1: K > 1 at Re=6400 t=40 ───────────────────────────────
    K_t40 = p['K']['Re6400_t40']
    certs.append(("K_t40_gt_1",
        f"K(Re=6400, t=40) = {K_t40} > 1",
        [f"(declare-const K Real)",
         f"(assert (= K {K_t40}))",
         f"(assert (> K 1.0))"]))

    # ── CERT-2: K crosses 1 (K at t≈5 < 1, K at t≈14 > 1) ──────────
    steps = p['staircase_growth']['steps']
    K_early = next(s['K_Lambda_eta'] for s in steps if abs(s['t'] - 5.22) < 1)
    K_cross = next(s['K_Lambda_eta'] for s in steps if abs(s['t'] - 14.0) < 1)
    certs.append(("K_crossover",
        f"K crosses 1: K(t≈5)={K_early} < 1 and K(t≈14)={K_cross} > 1",
        [f"(declare-const K_early Real)",
         f"(declare-const K_late Real)",
         f"(assert (= K_early {K_early}))",
         f"(assert (= K_late {K_cross}))",
         f"(assert (< K_early 1.0))",
         f"(assert (> K_late 1.0))"]))

    # ── CERT-3: Convergence ratio < 1.3 ─────────────────────────────
    r256 = p['convergence']['Re6400_t20']['ratio_256_128']
    certs.append(("convergence_ratio_bounded",
        f"ratio(256/128) = {r256} < 1.3",
        [f"(declare-const r Real)",
         f"(assert (= r {r256}))",
         f"(assert (< r 1.3))",
         f"(assert (> r 1.0))"]))

    # ── CERT-4: Convergence ratios decreasing (trending → 1) ────────
    r128 = p['convergence']['Re6400_t20']['ratio_128_64']
    certs.append(("convergence_ratios_decreasing",
        f"ratio(128/64)={r128} > ratio(256/128)={r256} (trending toward 1)",
        [f"(declare-const r1 Real)",
         f"(declare-const r2 Real)",
         f"(assert (= r1 {r128}))",
         f"(assert (= r2 {r256}))",
         f"(assert (> r1 r2))",
         f"(assert (> r2 1.0))"]))

    # ── CERT-5: Onset table final/trough is strictly increasing ──────
    onset = p['onset_table']
    decls = []
    asserts = []
    for i, row in enumerate(onset):
        decls.append(f"(declare-const ft_{i} Real)")
        asserts.append(f"(assert (= ft_{i} {row['final_over_trough']}))")
    for i in range(len(onset) - 1):
        asserts.append(f"(assert (< ft_{i} ft_{i+1}))")
    certs.append(("onset_ft_strictly_increasing",
        "final/trough strictly increases with Re across onset table",
        decls + asserts))

    # ── CERT-6: Staircase peaks increasing ───────────────────────────
    # Extract local peaks: t=20, t=30, t=38.8, t=39.92
    peak_times = [20.0, 30.0, 38.82, 39.92]
    peaks = []
    for target in peak_times:
        best = min(steps, key=lambda s: abs(s['t'] - target))
        peaks.append(best)
    decls = []
    asserts = []
    for i, pk in enumerate(peaks):
        decls.append(f"(declare-const peak_{i} Real)")
        asserts.append(f"(assert (= peak_{i} {pk['K_Lambda_eta']}))")
    for i in range(len(peaks) - 1):
        asserts.append(f"(assert (< peak_{i} peak_{i+1}))")
    certs.append(("staircase_peaks_increasing",
        "K at successive peaks increases (ratchet property)",
        decls + asserts))

    # ── CERT-7: Lambda grows while enstrophy decreases ───────────────
    L_t20 = p['Lambda']['Re6400_N128_t20']
    L_t40 = p['Lambda']['Re6400_N128_t40']
    O_t20 = p['enstrophy']['Re6400_t20']
    O_t40 = p['enstrophy']['Re6400_t40']
    certs.append(("Lambda_grows_Omega_decays",
        f"Λ: {L_t20}→{L_t40} (grows), Ω: {O_t20}→{O_t40} (decays)",
        [f"(declare-const L20 Real)", f"(declare-const L40 Real)",
         f"(declare-const O20 Real)", f"(declare-const O40 Real)",
         f"(assert (= L20 {L_t20}))", f"(assert (= L40 {L_t40}))",
         f"(assert (= O20 {O_t20}))", f"(assert (= O40 {O_t40}))",
         f"(assert (> L40 L20))",
         f"(assert (< O40 O20))"]))

    # ── CERT-8: Subcritical tameness ─────────────────────────────────
    sub = [r for r in onset if r['Re'] <= 800]
    decls = []
    asserts = []
    for i, row in enumerate(sub):
        decls.append(f"(declare-const ft_sub_{i} Real)")
        asserts.append(f"(assert (= ft_sub_{i} {row['final_over_trough']}))")
        asserts.append(f"(assert (<= ft_sub_{i} 1.1))")
    certs.append(("subcritical_tame",
        "Re ≤ 800 → final/trough ≤ 1.1",
        decls + asserts))

    # ── CERT-9: K definition consistency (spot check) ────────────────
    # K = Lambda * eta for Re=6400 t=40.  Pre-compute the product to stay in QF_LRA.
    K_computed = L_t40 * p['eta']['Re6400_t40']
    certs.append(("K_def_consistent",
        f"K = Λ·η: {L_t40} * {p['eta']['Re6400_t40']} = {K_computed:.4f} ≈ {K_t40}",
        [f"(declare-const K_computed Real)",
         f"(declare-const K_reported Real)",
         f"(assert (= K_computed {K_computed}))",
         f"(assert (= K_reported {K_t40}))",
         f"; Difference < 0.5 (rational approximation tolerance)",
         f"(assert (< (- K_computed K_reported) 0.5))",
         f"(assert (> (- K_computed K_reported) (- 0.5)))"]))

    # ── CERT-10: Regime boundaries consistent ────────────────────────
    certs.append(("regime_boundaries_ordered",
        "0 < regime1_upper < regime2_lower < regime2_upper < regime3_lower",
        [f"(declare-const r1u Real)", f"(declare-const r2l Real)",
         f"(declare-const r2u Real)", f"(declare-const r3l Real)",
         f"(assert (= r1u {p['Re_bifurcation']['regime_1_upper']}))",
         f"(assert (= r2l {p['Re_bifurcation']['regime_2_range'][0]}))",
         f"(assert (= r2u {p['Re_bifurcation']['regime_2_range'][1]}))",
         f"(assert (= r3l {p['Re_bifurcation']['regime_3_lower']}))",
         f"(assert (> r1u 0))",
         f"(assert (< r1u r2l))",
         f"(assert (< r2l r2u))",
         f"(assert (< r2u r3l))"]))

    # ── Write combined SMT2 file ─────────────────────────────────────
    for name, desc, smt_lines in certs:
        blocks.append(f"; === CERT: {name} ===")
        blocks.append(f"; {desc}")
        blocks.append("(push)")
        for line in smt_lines:
            blocks.append(line)
        blocks.append("(check-sat)")
        blocks.append("(pop)")
        blocks.append("")

    smt2_path = Z3_DIR / "all_certs.smt2"
    smt2_path.write_text("\n".join(blocks) + "\n")

    # ── Write runner script ──────────────────────────────────────────
    runner = Z3_DIR / "run_z3.py"
    cert_names = [c[0] for c in certs]
    runner.write_text(f'''#!/usr/bin/env python3
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

cert_names = {cert_names!r}

def run_z3_cli():
    r = subprocess.run(["z3", str(smt2)], capture_output=True, text=True)
    lines = [l.strip() for l in r.stdout.strip().split("\\n") if l.strip()]
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
        s.from_string("(set-logic QF_LRA)\\n" + inner)
        r = s.check()
        results.append("sat" if r == z3.sat else "unsat" if r == z3.unsat else "unknown")
    return results

if HAS_Z3_PY:
    lines = run_z3_python()
else:
    lines = run_z3_cli()

results = {{}}
for i, name in enumerate(cert_names):
    status = lines[i] if i < len(lines) else "error"
    results[name] = {{"status": status, "index": i}}
    mark = "✓" if status == "sat" else "✗"
    print(f"  {{mark}} {{name}}: {{status}}")

out = Z3_DIR / "results.json"
out.write_text(json.dumps(results, indent=2))
print(f"\\nWrote {{out}}")

failed = [n for n, r in results.items() if r["status"] != "sat"]
if failed:
    print(f"\\nFAILED: {{failed}}")
    sys.exit(1)
else:
    print(f"\\nAll {{len(results)}} certificates SAT")
''')
    runner.chmod(0o755)

    print(f"Wrote {smt2_path}  ({len(certs)} certificates)")
    print(f"Wrote {runner}")
    for name, desc, _ in certs:
        print(f"  CERT: {name}")


if __name__ == "__main__":
    main()
