#!/usr/bin/env python3
"""Single-snapshot crease visualization at N=256, t=28.

Renders the direction-field crease at the highest available resolution:
  - Volume render of omega_mag
  - Isosurface of T2_norm at 50% max (crease geometry)
  - Slice at z=pi with xi direction glyphs

Usage:
  pvbatch view_crease_t28.py
  paraview --script=view_crease_t28.py

Environment variables:
  SNAPSHOT_DIR   Path to HDF5/XDMF directory. Default: ../data/snapshots
  OUTPUT_DIR     Where to save PNG. Default: output/
"""

import os
import sys
import math

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SNAPSHOT_DIR = os.environ.get("SNAPSHOT_DIR",
    os.path.join(SCRIPT_DIR, "..", "data", "snapshots"))
OUTPUT_DIR = os.environ.get("OUTPUT_DIR",
    os.path.join(SCRIPT_DIR, "output"))

XDMF_FILE = os.path.join(SNAPSHOT_DIR, "run_256_t0028.0.xdmf")

def build_scene():
    try:
        from paraview.simple import (
            OpenDataFile, GetActiveView, CreateRenderView,
            Show, Hide, Threshold, Slice, Glyph,
            GetColorTransferFunction, GetOpacityTransferFunction,
            ColorBy, Text, SaveScreenshot, ResetCamera,
        )
    except ImportError:
        print("ERROR: Run inside ParaView (pvbatch or paraview --script)")
        sys.exit(1)

    if not os.path.exists(XDMF_FILE):
        print(f"ERROR: {XDMF_FILE} not found")
        print("Place HDF5 snapshots in data/snapshots/ or set SNAPSHOT_DIR")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Use existing view in interactive mode, create new only in pvbatch
    view = GetActiveView()
    if view is None:
        view = CreateRenderView()
    view.ViewSize = [1200, 1000]
    view.Background = [0.05, 0.05, 0.1]

    reader = OpenDataFile(XDMF_FILE)
    reader.UpdatePipeline()

    # Volume render omega_mag
    disp = Show(reader, view)
    disp.Representation = "Volume"
    ColorBy(disp, ("POINTS", "omega_mag"))
    lut = GetColorTransferFunction("omega_mag")
    lut.ApplyPreset("Viridis (matplotlib)", True)

    pwf = GetOpacityTransferFunction("omega_mag")
    data_range = reader.PointData["omega_mag"].GetRange()
    span = data_range[1] - data_range[0]
    if span > 0:
        pwf.Points = [
            data_range[0], 0.0, 0.5, 0.0,
            data_range[0] + 0.25 * span, 0.01, 0.5, 0.0,
            data_range[0] + 0.5 * span, 0.05, 0.5, 0.0,
            data_range[1], 0.25, 0.5, 0.0,
        ]

    # T2_norm crease isosurface
    try:
        t2_range = reader.PointData["T2_norm"].GetRange()
        thresh = Threshold(Input=reader)
        thresh.Scalars = ["POINTS", "T2_norm"]
        thresh.LowerThreshold = 0.5 * t2_range[1]
        thresh.UpperThreshold = t2_range[1]
        thresh.UpdatePipeline()
        td = Show(thresh, view)
        td.Representation = "Surface"
        td.DiffuseColor = [1.0, 0.2, 0.05]
        td.Opacity = 0.35
    except Exception as e:
        print(f"T2_norm threshold skipped: {e}")

    # Slice + xi glyphs
    try:
        sl = Slice(Input=reader)
        sl.SliceType = "Plane"
        sl.SliceType.Origin = [math.pi, math.pi, math.pi]
        sl.SliceType.Normal = [0, 0, 1]
        sl.UpdatePipeline()

        gl = Glyph(Input=sl)
        gl.GlyphType = "Arrow"
        gl.OrientationArray = ["POINTS", "xi"]
        gl.ScaleArray = ["POINTS", "No scale array"]
        gl.ScaleFactor = 2.0 * math.pi / 256 * 8
        gl.MaximumNumberOfSamplePoints = (256 // 8) ** 2
        gl.GlyphMode = "Every Nth Point"
        gl.Stride = 8
        gl.UpdatePipeline()

        gd = Show(gl, view)
        gd.Representation = "Surface"
        gd.DiffuseColor = [0.9, 0.9, 0.9]
        gd.Opacity = 0.6
        Hide(sl, view)
    except Exception as e:
        print(f"xi glyph skipped: {e}")

    # Annotation
    txt = Text()
    txt.Text = "N=256  t=28.0  Re=6400"
    td = Show(txt, view)
    td.FontSize = 20
    td.Color = [1, 1, 1]
    td.WindowLocation = "Upper Center"

    ResetCamera(view)
    from paraview.simple import Render
    Render(view)

    # Save PNG only in pvbatch (non-interactive) mode
    interactive = True
    try:
        from paraview.simple import servermanager
        interactive = servermanager.vtkProcessModule.GetProcessModule().GetPartitionId() == 0 \
            and not servermanager.vtkProcessModule.GetProcessModule().GetOptions().GetSymmetricMPIMode()
        # More reliable check: pvbatch sets sys.argv[0] to pvbatch path
        interactive = "pvbatch" not in sys.argv[0].lower()
    except Exception:
        pass

    if not interactive or "--save" in sys.argv:
        out = os.path.join(OUTPUT_DIR, "crease_256_t28.png")
        SaveScreenshot(out, view=view, ImageResolution=[1200, 1000])
        print(f"Saved: {out}")
    else:
        print("Scene loaded. Interact in the GUI. Use --save or pvbatch to export PNG.")


if __name__ == "__main__":
    print(f"view_crease_t28.py")
    print(f"XDMF: {XDMF_FILE}")
    build_scene()
