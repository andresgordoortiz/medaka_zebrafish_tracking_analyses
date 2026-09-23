"""Surface flow maps for medaka 3005 at the equivalent key timepoints.

Adapted from mk2508_surface_flow_maps.py for the 3005 dataset structure:
  - oriented_tracks_medaka.csv (single per-spot CSV, columns: TRACK_ID,
    FRAME, POSITION_X, POSITION_Y, POSITION_Z, RADIAL_DIST,
    SPHERICAL_DEPTH, THETA_DEG, PHI_DEG)
  - sphere_params.csv
  - tracks.csv (the raw ultrack output, only TRACK_ID column with index)

The 3005 dataset is the HIGH-RESOLUTION medaka acquisition (different
voxel size than the low-res 2508 dataset) so we override the voxel and
the surface/cross-section parameters below.

Voxel calibration (3005, high-res)
---------------------------------
  voxel_size = {"x_um": 0.5272727272727272,
                "y_um": 0.5272727272727272,
                "z_um": 0.5263157894736842}
  shape    = (T=395, Z=266, Y=760, X=760)
  frame interval = 30 s (assumed; check dataset metadata if available)

Why a TAILORED cross-section
----------------------------
The 3005 sphere fit gives R = 803.99 px = 423.9 µm (vs 840.0 µm for 2508).
The 3005 embryo is physically about half the size of the 2508 embryo, so
we scale every spatial parameter by ~1/2 to keep the figure framing
comparable:

  2508 (low-res)  | 3005 (high-res) | rationale
  --------------- + --------------- + ------------------------------
  VOXEL: 1.05152   | VOXEL: 0.52727  | per metadata
  SURFACE: 30 µm   | SURFACE: 15 µm  | half-thickness surface band
  FLOW_BIN: 30 µm  | FLOW_BIN: 15 µm | half-size spatial bin
  Z_BAND: 500 µm   | Z_BAND: 250 µm  | half-FOV projection depth

Frame selection
---------------
Anchors (user-confirmed biological landmarks):
  - Start of internalisation: frame 180
  - First contraction (end):  frame 395
At 30 s/frame (2 frames/min):
  35   - 3 h before first contraction (395 - 360 frames)
  155  - 2 h before first contraction (395 - 240 frames)
  180  - Start of internalisation
  290  - Mid-epiboly / pre-contraction (midway between 180 and 395)
  395  - First contraction (end)

Run with:
  uv run python results/medaka_dynamics/mk3005_surface_flow_maps.py
"""

from __future__ import annotations

# Reuse the analysis logic from the 2508 script: import it directly so any
# fixes propagate to both datasets.
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))

# Import the heavy lifting from the 2508 script.
import mk2508_surface_flow_maps as m
from mk2508_surface_flow_maps import (
    apply_voxel_calibration,
    compute_per_step_velocity,
    compute_track_stats_per_frame,
    select_tracks,
    select_window,
    select_surface,
    select_slice,
    bin_xy,
    classify_bins,
    draw_panel,
    _build_page1,
    _build_page2,
    FRAME_INTERVAL_SEC,
    FRAMES_PER_MIN,
    MAX_INSTANT_SPEED_UM_PER_FRAME,
    MIN_TRACK_LEN,
    ARROW_HEAD_FRAC,
    ARROW_LW,
    Z_MID,
    EXTENT_PAD_UM,
    FRAME_LABELS,
)

# ---- 3005-specific overrides (high-resolution acquisition) -------------
# Voxel (from the metadata the user provided)
m.VOXEL_XY_UM_PX = 0.5272727272727272    # µm per XY pixel
m.VOXEL_Z_UM_PX  = 0.5263157894736842    # µm per Z pixel
# Cross-section tailored for the smaller 3005 embryo (R ≈ 424 µm vs 840 µm
# for 2508; see module docstring above).
m.SURFACE_UM = 15.0    # spherical depth <= this counts as "surface"
m.FLOW_BIN_UM = 15.0   # spatial bin in the cross-section
m.Z_BAND_UM  = 250.0   # |z - Z_MID| <= this band (full FOV projection)
m.FLOW_MIN_N = 4       # minimum # vectors per bin to draw an arrow
m.HALF_WIN   = 8       # frames before/after for velocity estimation
# Arrow scale: arrow length per µm/frame.  Same numerical value as 2508
# (50) but represents half the absolute displacement, so the visual
# density matches the 2508 panel.
m.ARROW_SCALE = 50
# Re-bind the local aliases so the rest of the script reads consistently.
VOXEL_XY_UM_PX = m.VOXEL_XY_UM_PX
VOXEL_Z_UM_PX  = m.VOXEL_Z_UM_PX
SURFACE_UM     = m.SURFACE_UM
FLOW_BIN_UM    = m.FLOW_BIN_UM
Z_BAND_UM      = m.Z_BAND_UM
FLOW_MIN_N     = m.FLOW_MIN_N
HALF_WIN       = m.HALF_WIN
ARROW_SCALE    = m.ARROW_SCALE

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# ---- Key frames for the 3005 dataset.
KEY_FRAMES_3005 = [35, 155, 180, 290, 395]
FRAME_LABELS_3005 = {
    35:  "3 h before first contraction",
    155: "2 h before first contraction",
    180: "Start of internalisation",
    290: "Mid-epiboly / pre-contraction",
    395: "First contraction (end)",
}


def load_3005(io_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return (oriented_spots in px, per-track-length table, sphere_params).

    The 3005 dataset doesn't ship a ``filtered_tracks.csv`` so we synthesise
    a per-track-length table directly from the oriented tracks file.  The
    track-level speed filter would normally come from the upstream file;
    here we just use the track length as the only track-level gate (which
    matches the 2508 pipeline because the data-derived 3-sigma speed
    threshold on filtered_tracks dropped zero tracks anyway).
    """
    spots = pd.read_csv(io_dir / "oriented_tracks_medaka.csv", low_memory=False)
    sph = pd.read_csv(io_dir / "sphere_params.csv")

    # Build a track-length table from the oriented spots.
    tl = (
        spots.groupby("TRACK_ID", as_index=False)
        .size()
        .rename(columns={"size": "NUMBER_SPOTS"})
    )

    numeric_keys = {"radius", "center_x", "center_y", "center_z",
                    "rmse", "cap_height", "cap_base_radius"}
    sphere = {row["parameter"]: row["value"] for _, row in sph.iterrows()}
    for k in list(sphere.keys()):
        if k in numeric_keys:
            try:
                sphere[k] = float(sphere[k])
            except (TypeError, ValueError):
                pass
    return spots, tl, sphere


def select_tracks_3005(track_table: pd.DataFrame,
                       min_len: int = MIN_TRACK_LEN) -> tuple[pd.Series, dict]:
    """Select tracks with at least ``min_len`` spots (no speed filter — see note)."""
    n0 = track_table["TRACK_ID"].nunique()
    long_enough = track_table["NUMBER_SPOTS"] >= min_len
    keep = track_table[long_enough]["TRACK_ID"]
    diag = {
        "n_total": n0,
        "n_long": int(long_enough.sum()),
        "n_kept": int(len(keep)),
    }
    return keep, diag


def main() -> None:
    # 3005 data lives in data/oriented_medaka_ultrack/, not in results/medaka_dynamics/
    io_dir = REPO / "data/oriented_medaka_ultrack"
    save_stem = HERE / "mk3005_surface_flow_maps"
    frames = KEY_FRAMES_3005

    print(f"\n{'='*70}\n  MEDAKA 3005 SURFACE FLOW MAPS\n{'='*70}\n")

    print("Loading data...")
    spots_px, track_table, sphere_px = load_3005(io_dir)
    print(f"  oriented_tracks_medaka : {len(spots_px):,} rows (raw, in pixels)")
    print(f"  sphere center          : ({sphere_px['center_x']:.1f}, "
          f"{sphere_px['center_y']:.1f}, {sphere_px['center_z']:.1f}) px  "
          f"R = {sphere_px['radius']:.1f} px")

    # Voxel calibration
    print(f"\nVoxel calibration: {VOXEL_XY_UM_PX} µm/px (XY), {VOXEL_Z_UM_PX} µm/px (Z), "
          f"{FRAME_INTERVAL_SEC}s/frame -> {FRAMES_PER_MIN} frames/min")
    spots, sphere = apply_voxel_calibration(spots_px, sphere_px)
    print(f"  sphere R in µm: {sphere['radius']:.1f} µm")
    print(f"  cell range X (µm): [{spots['POSITION_X'].min():.1f}, {spots['POSITION_X'].max():.1f}]")
    print(f"  cell range Y (µm): [{spots['POSITION_Y'].min():.1f}, {spots['POSITION_Y'].max():.1f}]")

    # Track-length filter
    keep, diag = select_tracks_3005(track_table)
    print("\nTrack filtering:")
    print(f"  total tracks         : {diag['n_total']:,}")
    print(f"  with n_spots >= {MIN_TRACK_LEN}  : {diag['n_long']:,}")
    print(f"  tracks kept          : {diag['n_kept']:,}  "
          f"({100 * diag['n_kept'] / diag['n_total']:.1f}%)")

    # Per-step velocity
    print("\nComputing per-step velocity (one row per step in a track)...")
    spots_v_all = compute_per_step_velocity(spots)
    n_dropped_speed = spots_v_all.attrs.get("n_dropped_speed", 0)
    spots_v = spots_v_all[spots_v_all["TRACK_ID"].isin(keep)]
    print(f"  per-step rows (after kept-track filter)     : {len(spots_v):,}")
    print(f"  per-step rows dropped (instant speed > "
          f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.1f} µm/frame "
          f"= {MAX_INSTANT_SPEED_UM_PER_FRAME * m.VOXEL_XY_UM_PX:.2f} µm/frame): "
          f"{n_dropped_speed:,} "
          f"({100 * n_dropped_speed / max(len(spots_v_all), 1):.2f}%)")
    print(f"  instantaneous XY speed stats (µm / min, after voxel + per-step filter):")
    print(f"    median = {spots_v['speed'].median() * FRAMES_PER_MIN:.3f}   "
          f"mean = {spots_v['speed'].mean() * FRAMES_PER_MIN:.3f}   "
          f"std = {spots_v['speed'].std() * FRAMES_PER_MIN:.3f}")

    # Build page 1 (flow maps) -- we have to patch FRAME_LABELS for the 3005
    # set so the panel titles show the right biological label.
    from mk2508_surface_flow_maps import FRAME_LABELS as LABELS_2508
    LABELS_2508.clear()
    LABELS_2508.update(FRAME_LABELS_3005)

    print(f"\nGenerating flow maps for frames {frames} ...")

    # Manually iterate the same way as make_figure() in the 2508 module so we
    # get the per-frame summaries; this avoids re-importing FRAME_LABELS.
    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 6.5))
    if n == 1:
        axes = [axes]
    summaries = []
    for ax, frame in zip(axes, frames):
        window = select_window(spots_v, frame, keep, HALF_WIN)
        surface = select_surface(window, SURFACE_UM)
        slice_ = select_slice(surface, Z_MID, Z_BAND_UM)
        fl = bin_xy(slice_, FLOW_BIN_UM, FLOW_MIN_N)
        fl = classify_bins(fl)
        bg = slice_[["MID_X", "MID_Y"]].to_numpy()
        draw_panel(ax, fl, sphere,
                   frame=frame, surface_um=SURFACE_UM,
                   z_mid=Z_MID, z_band=Z_BAND_UM,
                   bin_um=FLOW_BIN_UM, arrow_scale=ARROW_SCALE,
                   arrow_lw=ARROW_LW, arrow_head=ARROW_HEAD_FRAC,
                   bg_scatter=bg)
        if len(fl) > 0:
            v_mean = fl["speed"].mean()
            up_frac = (fl["move_type"] == "AP-ward").mean()
            vx_mean = fl["vx"].mean()
            vy_mean = fl["vy"].mean()
        else:
            v_mean = vx_mean = vy_mean = up_frac = np.nan
        summaries.append({
            "frame": frame,
            "label": FRAME_LABELS_3005.get(frame, ""),
            "n_vectors": int(len(slice_)),
            "n_bins": int(len(fl)),
            "n_surface_window": int(len(surface)),
            "mean_speed_um_per_frame": float(v_mean),
            "mean_speed_um_per_min":   float(v_mean * FRAMES_PER_MIN),
            "mean_vx_um_per_frame": float(vx_mean),
            "mean_vy_um_per_frame": float(vy_mean),
            "frac_ap_ward_bins": float(up_frac),
        })
    handles = [
        Line2D([0], [0], color="#2166AC", lw=2,
               label="AP-ward (vy < 0, toward animal pole)"),
        Line2D([0], [0], color="#B2182B", lw=2,
               label="VP-ward (vy > 0, toward vegetal pole)"),
        Line2D([0], [0], color="black", lw=1,
               label=f"Sphere (R = {sphere['radius']:.0f} µm)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), frameon=False)
    speed_clause = (
        f"tracks ≥ {MIN_TRACK_LEN} spots; per-step instant speed ≤ "
        f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.1f} µm/frame "
        f"({n_dropped_speed:,} comet steps dropped)"
    )
    fig.suptitle(
        "Medaka 3005 — Surface flow maps at equivalent key timepoints\n"
        f"Side view (XY, AP at top) | surface band depth ≤ {SURFACE_UM:.0f} µm | "
        f"voxel {VOXEL_XY_UM_PX} µm/px | {speed_clause}",
        fontsize=11, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    page1_fig = fig
    plt.close(fig)

    # Print summary table
    summary_df = pd.DataFrame(summaries)
    print("\nSummary table:")
    print(summary_df.to_string(index=False))

    # Save summary CSV
    summary_path = io_dir / "mk3005_surface_flow_maps_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    # Build the 2-page PDF — Panel B from BIN-AVERAGED stats (matches figure)
    pdf_path = save_stem.with_suffix(".pdf")
    track_stats = compute_track_stats_per_frame(
        spots_v=spots_v, valid_tracks=keep, frames=frames,
        half_win=HALF_WIN, surface_um=SURFACE_UM,
        bin_um=FLOW_BIN_UM, min_n=FLOW_MIN_N,
    )
    print("\nBin-averaged Panel B statistics (matches the figure):")
    for f in frames:
        st = track_stats[f]
        print(f"  frame {f}: "
              f"frac AP-ward bins = {st['frac_ap_ward_bins']*100:.1f}%   "
              f"mean bin speed = {st['mean_bin_speed_um_per_min']:.3f} µm/min   "
              f"({st['n_bins']:,} bins, {st['n_steps']:,} steps)")

    with PdfPages(pdf_path) as pdf:
        pdf.savefig(page1_fig, bbox_inches="tight")
        plt.close(page1_fig)
        page2 = _build_page2(summaries, track_summary=track_stats)
        pdf.savefig(page2, bbox_inches="tight")
        plt.close(page2)

    # Save page 1 PNG separately
    page1_fig2 = _build_page1(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=frames, surface_um=SURFACE_UM, z_mid=Z_MID, z_band=Z_BAND_UM,
        bin_um=FLOW_BIN_UM, min_n=FLOW_MIN_N,
        arrow_scale=ARROW_SCALE, arrow_lw=ARROW_LW,
        arrow_head=ARROW_HEAD_FRAC, half_win=HALF_WIN,
        speed_thr_px=None, speed_median_px=None,
        speed_std_px=None, n_dropped_speed=n_dropped_speed,
    )
    page1_fig2.savefig(save_stem.with_suffix(".png"), dpi=130, bbox_inches="tight")
    plt.close(page1_fig2)

    print(f"\nSaved 2-page PDF to {pdf_path}")
    print(f"Saved page-1 PNG to {save_stem}.png")
    print(f"Saved summary to {summary_path}")


if __name__ == "__main__":
    main()