#!/usr/bin/env python3
"""S6: Extract proof gaps from lake build errors.

Parses `lake build` stderr, identifies sorry obligations and type errors,
and generates standalone lemma stubs in generated/Gaps.lean.

Usage: lake build 2>&1 | python3 scripts/s6_gap_extract.py
   or: python3 scripts/s6_gap_extract.py  (runs lake build internally)
"""

import re
import subprocess
import sys
from pathlib import Path

OUT = Path(__file__).parent.parent / "generated" / "Gaps.lean"


def run_build():
    """Run lake build and capture stderr."""
    r = subprocess.run(
        ["lake", "build"],
        capture_output=True, text=True,
        cwd=Path(__file__).parent.parent,
    )
    return r.stdout + "\n" + r.stderr


def parse_errors(output):
    """Extract sorry locations and unsolved goals."""
    gaps = []

    # Pattern: file:line:col: warning: declaration uses 'sorry'
    sorry_pat = re.compile(r"(.+?):(\d+):(\d+):\s*warning.*sorry")
    for m in sorry_pat.finditer(output):
        gaps.append({
            "file": m.group(1),
            "line": int(m.group(2)),
            "type": "sorry",
            "msg": m.group(0).strip(),
        })

    # Pattern: unsolved goals
    goal_pat = re.compile(r"(.+?):(\d+):(\d+):\s*error.*unsolved goals?\s*\n((?:.*\n)*?)(?=\S+:\d+:|$)")
    for m in goal_pat.finditer(output):
        gaps.append({
            "file": m.group(1),
            "line": int(m.group(2)),
            "type": "unsolved_goal",
            "msg": m.group(4).strip(),
        })

    # Pattern: type mismatch
    type_pat = re.compile(r"(.+?):(\d+):(\d+):\s*error.*type mismatch")
    for m in type_pat.finditer(output):
        gaps.append({
            "file": m.group(1),
            "line": int(m.group(2)),
            "type": "type_mismatch",
            "msg": m.group(0).strip(),
        })

    return gaps


def generate_gaps_lean(gaps):
    """Write standalone lemma stubs."""
    lines = []
    lines.append("/-!")
    lines.append("  Gaps.lean — Auto-generated from lake build errors")
    lines.append(f"  {len(gaps)} gaps found")
    lines.append("-/")
    lines.append("")

    sorry_gaps = [g for g in gaps if g["type"] == "sorry"]
    error_gaps = [g for g in gaps if g["type"] != "sorry"]

    lines.append(f"-- {len(sorry_gaps)} sorry obligations")
    lines.append(f"-- {len(error_gaps)} type/goal errors")
    lines.append("")

    for i, g in enumerate(gaps):
        lines.append(f"-- GAP {i+1}: {g['type']} at {g['file']}:{g['line']}")
        lines.append(f"-- {g['msg'][:200]}")
        lines.append("")

    return "\n".join(lines) + "\n"


def main():
    if not sys.stdin.isatty():
        output = sys.stdin.read()
    else:
        print("Running lake build...")
        output = run_build()

    gaps = parse_errors(output)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(generate_gaps_lean(gaps))

    print(f"Wrote {OUT}")
    print(f"  {len(gaps)} gaps total")
    for g in gaps:
        print(f"  [{g['type']}] {g['file']}:{g['line']}")

    # Also print sorry count
    sorry_count = sum(1 for g in gaps if g["type"] == "sorry")
    print(f"\nSorry count: {sorry_count}")


if __name__ == "__main__":
    main()
