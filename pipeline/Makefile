.PHONY: install params s1 s7 z3-run s8 build gaps loop status clean all

# ── Prerequisites ────────────────────────────────────────────────────
install:
	@echo "=== Installing elan (Lean 4) ==="
	curl https://elan-init.lean-lang.org/elan-init.sh -sSf | sh -s -- -y --default-toolchain leanprover/lean4:stable
	@echo "=== Installing Z3 (Python) ==="
	pip install z3-solver
	@echo "=== Verify ==="
	lean --version
	python3 -c "import z3; print('z3', z3.get_version_string())"

# ── S0: Extract params from V5 databases ─────────────────────────────
params:
	python3 scripts/extract_params.py

# ── S1: Generate Lean definitions ────────────────────────────────────
s1: params
	python3 scripts/s1_gen_lean_defs.py

# ── S7: Generate Z3 certificates ────────────────────────────────────
s7: params
	python3 scripts/s7_z3_certs.py

# ── Z3: Run certificates ────────────────────────────────────────────
z3-run: s7
	python3 z3/run_z3.py

# ── S8: Bridge Z3 → Lean axioms ─────────────────────────────────────
s8: z3-run
	python3 scripts/s8_z3_to_lean.py

# ── Build: lake build ────────────────────────────────────────────────
build:
	lake build 2>&1 | tee build.log
	@echo "=== Build complete ==="

# ── S6: Extract gaps ────────────────────────────────────────────────
gaps:
	lake build 2>&1 | python3 scripts/s6_gap_extract.py

# ── Loop: build + extract gaps ───────────────────────────────────────
loop: gaps
	@echo "Fix gaps in generated/Gaps.lean, then run 'make loop' again"

# ── Generate all (no build) ─────────────────────────────────────────
all: s1 s8
	@echo "=== All generated files ready ==="
	@echo "Generated:"
	@ls -la generated/
	@echo ""
	@echo "Next: configure lakefile.lean, then 'make build'"

# ── Status ───────────────────────────────────────────────────────────
status:
	@echo "=== sorry count ==="
	@grep -c sorry NSLean4/NSTheorems.lean 2>/dev/null || echo "0"
	@echo "=== Z3 status ==="
	@python3 -c "import json; r=json.load(open('z3/results.json')); sat=sum(1 for v in r.values() if v['status']=='sat'); print(f'{sat}/{len(r)} SAT')" 2>/dev/null || echo "not run"
	@echo "=== generated files ==="
	@ls -la generated/ 2>/dev/null || echo "none"

# ── Clean ────────────────────────────────────────────────────────────
clean:
	rm -rf generated/ z3/results.json build.log .lake/build/
