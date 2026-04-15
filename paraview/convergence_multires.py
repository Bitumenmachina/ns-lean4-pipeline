#!/usr/bin/env python3
"""Multi-resolution convergence visualization for ParaView.

Loads N=64, N=128, N=256 XDMF snapshots at a matched timestep and renders
three side-by-side viewports showing:
  - Volume render of omega_mag (viridis)
  - Threshold isosurface on T2_norm at 50% of per-resolution max (crease)
  - Slice at z=pi with xi direction glyphs (stride 8)

Annotation overlay per viewport: "N=XX, Lambda_max=YYY, K=ZZ"
Values pulled from convergence_summary.csv at matched t.

Usage:
  # Headless (saves PNG):
  pvbatch convergence_multires.py [--time 20]

  # Interactive:
  paraview --script=convergence_multires.py

Environment variables:
  SNAPSHOT_DIR   Path to directory containing run_*_tXXXX.xdmf files
                 Default: ../data/snapshots (relative to this script)
  CONVERGENCE_CSV  Path to convergence_summary.csv
                 Default: ../data/convergence/convergence_summary.csv
  OUTPUT_DIR     Where to save PNGs. Default: output/ (relative to this script)
  TIME           Timestep to visualize. Default: 20. Overridden by --time arg.

Requirements:
  - ParaView >= 5.10 with Python scripting
  - XDMF files + companion HDF5 files in SNAPSHOT_DIR
  - convergence_summary.csv in CONVERGENCE_CSV path
"""

import os
import sys
import csv
import math

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def parse_args():
    """Parse --time argument from sys.argv."""
    t = float(os.environ.get("TIME", "20"))
    for i, arg in enumerate(sys.argv):
        if arg == "--time" and i + 1 < len(sys.argv):
            t = float(sys.argv[i + 1])
    return t

TARGET_TIME = parse_args()
SNAPSHOT_DIR = os.environ.get("SNAPSHOT_DIR",
    os.path.join(SCRIPT_DIR, "..", "data", "snapshots"))
CONVERGENCE_CSV = os.environ.get("CONVERGENCE_CSV",
    os.path.join(SCRIPT_DIR, "..", "data", "convergence", "convergence_summary.csv"))
OUTPUT_DIR = os.environ.get("OUTPUT_DIR",
    os.path.join(SCRIPT_DIR, "output"))

RESOLUTIONS = [64, 128, 256]
VIEW_SIZE = (2400, 800)  # 3 panels at 800x800

# Viridis-inspired 4-point opacity ramp for volume rendering
OPACITY_POINTS = [
    (0.0, 0.0),
    (0.25, 0.02),
    (0.5, 0.08),
    (1.0, 0.3),
]

GLYPH_STRIDE = 8
SLICE_Z = math.pi  # z = pi midplane

# ---------------------------------------------------------------------------
# Load annotation data from CSV
# ---------------------------------------------------------------------------

def load_annotations(csv_path, target_t):
    """Load Lambda and K values from convergence_summary.csv at target_t."""
    annotations = {}
    if not os.path.exists(csv_path):
        print(f"WARNING: {csv_path} not found, annotations will be empty")
        return annotations

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t = float(row["t"])
            except (ValueError, KeyError):
                continue
            if abs(t - target_t) < 0.5:
                for n in RESOLUTIONS:
                    lam_key = f"\u039b_{n}"  # Λ_XX
                    k_key = f"K_{n}"
                    # Try both unicode and ASCII column names
                    lam_val = row.get(lam_key, row.get(f"Λ_{n}", "—"))
                    k_val = row.get(k_key, "—")
                    annotations[n] = {
                        "Lambda": lam_val,
                        "K": k_val,
                    }
                break
    return annotations


def find_xdmf(snapshot_dir, n, target_t):
    """Find the XDMF file closest to target_t for resolution n."""
    import glob
    pattern = os.path.join(snapshot_dir, f"run_{n}_t*.xdmf")
    candidates = glob.glob(pattern)
    if not candidates:
        return None

    def extract_time(path):
        base = os.path.basename(path)
        # run_64_t0020.0.xdmf -> 20.0
        t_str = base.split("_t")[1].replace(".xdmf", "")
        return float(t_str)

    best = min(candidates, key=lambda p: abs(extract_time(p) - target_t))
    return best


# ---------------------------------------------------------------------------
# ParaView rendering
# ---------------------------------------------------------------------------

def build_scene():
    """Build the three-viewport ParaView scene."""
    try:
        from paraview.simple import (
            OpenDataFile, GetActiveView, CreateRenderView,
            GetDisplayProperties, Show, Hide,
            Threshold, Slice, Glyph,
            GetColorTransferFunction, GetOpacityTransferFunction,
            GetLayout, AssignViewToLayout,
            Text, SaveScreenshot, ResetCamera,
            SetActiveView, GetActiveSource,
            ColorBy, servermanager,
        )
    except ImportError:
        print("ERROR: This script must be run inside ParaView (pvbatch or paraview --script)")
        print("Usage: pvbatch convergence_multires.py [--time 20]")
        sys.exit(1)

    annotations = load_annotations(CONVERGENCE_CSV, TARGET_TIME)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # In pvbatch, no layout exists yet — create the first view explicitly
    first_view = CreateRenderView()
    layout = GetLayout(first_view)
    views = []
    sources = []

    for idx, n in enumerate(RESOLUTIONS):
        xdmf_path = find_xdmf(SNAPSHOT_DIR, n, TARGET_TIME)
        if xdmf_path is None:
            print(f"WARNING: No XDMF found for N={n} at t={TARGET_TIME}")
            continue

        print(f"Loading N={n}: {xdmf_path}")

        # Create render view
        if idx == 0:
            view = first_view
        else:
            view = CreateRenderView()

        view.ViewSize = [VIEW_SIZE[0] // 3, VIEW_SIZE[1]]
        view.Background = [0.1, 0.1, 0.15]
        views.append(view)
        SetActiveView(view)

        # Load XDMF
        reader = OpenDataFile(xdmf_path)
        reader.UpdatePipeline()

        # --- Layer 1: Volume render of omega_mag ---
        disp = Show(reader, view)
        disp.Representation = "Volume"
        ColorBy(disp, ("POINTS", "omega_mag"))

        lut = GetColorTransferFunction("omega_mag")
        lut.ApplyPreset("Viridis (matplotlib)", True)

        pwf = GetOpacityTransferFunction("omega_mag")
        data_range = reader.PointData["omega_mag"].GetRange()
        if data_range[1] > data_range[0]:
            span = data_range[1] - data_range[0]
            points = []
            for frac, opacity in OPACITY_POINTS:
                points.extend([data_range[0] + frac * span, opacity, 0.5, 0.0])
            pwf.Points = points

        # --- Layer 2: Threshold isosurface on T2_norm (crease) ---
        try:
            t2_range = reader.PointData["T2_norm"].GetRange()
            t2_threshold = 0.5 * t2_range[1]
            if t2_threshold > 0:
                thresh = Threshold(Input=reader)
                thresh.Scalars = ["POINTS", "T2_norm"]
                thresh.LowerThreshold = t2_threshold
                thresh.UpperThreshold = t2_range[1]
                thresh.UpdatePipeline()

                thresh_disp = Show(thresh, view)
                thresh_disp.Representation = "Surface"
                thresh_disp.DiffuseColor = [1.0, 0.3, 0.1]
                thresh_disp.Opacity = 0.4
        except Exception as e:
            print(f"  T2_norm threshold skipped for N={n}: {e}")

        # --- Layer 3: Slice at z=pi with xi glyphs ---
        try:
            sl = Slice(Input=reader)
            sl.SliceType = "Plane"
            sl.SliceType.Origin = [math.pi, math.pi, SLICE_Z]
            sl.SliceType.Normal = [0, 0, 1]
            sl.UpdatePipeline()

            gl = Glyph(Input=sl)
            gl.GlyphType = "Arrow"
            gl.OrientationArray = ["POINTS", "xi"]
            gl.ScaleArray = ["POINTS", "No scale array"]
            gl.ScaleFactor = 2.0 * math.pi / n * GLYPH_STRIDE
            gl.MaximumNumberOfSamplePoints = (n // GLYPH_STRIDE) ** 2
            gl.GlyphMode = "Every Nth Point"
            gl.Stride = GLYPH_STRIDE
            gl.UpdatePipeline()

            gl_disp = Show(gl, view)
            gl_disp.Representation = "Surface"
            gl_disp.DiffuseColor = [0.9, 0.9, 0.9]
            gl_disp.Opacity = 0.7

            # Hide the slice itself
            Hide(sl, view)
        except Exception as e:
            print(f"  xi glyph skipped for N={n}: {e}")

        # --- Annotation ---
        ann = annotations.get(n, {"Lambda": "—", "K": "—"})
        txt = Text()
        txt.Text = f"N={n}  \\Lambda={ann['Lambda']}  K={ann['K']}"
        txt_disp = Show(txt, view)
        txt_disp.FontSize = 18
        txt_disp.Color = [1.0, 1.0, 1.0]
        txt_disp.WindowLocation = "Upper Center"

        ResetCamera(view)
        sources.append(reader)

    # --- Link cameras ---
    if len(views) >= 2:
        for v in views[1:]:
            v.CameraPosition = views[0].CameraPosition
            v.CameraFocalPoint = views[0].CameraFocalPoint
            v.CameraViewUp = views[0].CameraViewUp
            v.CameraViewAngle = views[0].CameraViewAngle

    # --- Save per-view images then composite ---
    panel_paths = []
    panel_w = VIEW_SIZE[0] // 3
    panel_h = VIEW_SIZE[1]

    for i, v in enumerate(views):
        n = RESOLUTIONS[i]
        panel_path = os.path.join(OUTPUT_DIR,
            f"_panel_N{n}_t{int(TARGET_TIME)}.png")
        SetActiveView(v)
        v.ViewSize = [panel_w, panel_h]
        try:
            SaveScreenshot(panel_path, view=v,
                ImageResolution=[panel_w, panel_h])
            panel_paths.append(panel_path)
            print(f"  Panel N={n}: {panel_path}")
        except Exception as e:
            print(f"  Could not save panel N={n}: {e}")

    # Composite panels side-by-side using PIL/Pillow or raw concat
    out_path = os.path.join(OUTPUT_DIR,
        f"convergence_t{int(TARGET_TIME)}.png")

    if len(panel_paths) >= 2:
        try:
            from PIL import Image
            panels = [Image.open(p) for p in panel_paths]
            total_w = sum(p.width for p in panels)
            max_h = max(p.height for p in panels)
            composite = Image.new("RGB", (total_w, max_h))
            x_offset = 0
            for p in panels:
                composite.paste(p, (x_offset, 0))
                x_offset += p.width
            composite.save(out_path)
            print(f"Saved composite: {out_path}")
            # Clean up panel files
            for p in panel_paths:
                os.remove(p)
        except ImportError:
            # No PIL — keep individual panels
            print(f"PIL not available; panels saved individually in {OUTPUT_DIR}/")
            # Rename first panel as the main output
            import shutil
            shutil.copy(panel_paths[0], out_path)
    elif len(panel_paths) == 1:
        import shutil
        shutil.copy(panel_paths[0], out_path)
        print(f"Saved (single panel): {out_path}")

    print(f"\nDone. Target time: t={TARGET_TIME}")
    print(f"Resolutions loaded: {[n for n in RESOLUTIONS if find_xdmf(SNAPSHOT_DIR, n, TARGET_TIME)]}")


if __name__ == "__main__":
    print(f"convergence_multires.py — t={TARGET_TIME}")
    print(f"SNAPSHOT_DIR: {SNAPSHOT_DIR}")
    print(f"CONVERGENCE_CSV: {CONVERGENCE_CSV}")
    print(f"OUTPUT_DIR: {OUTPUT_DIR}")
    build_scene()
