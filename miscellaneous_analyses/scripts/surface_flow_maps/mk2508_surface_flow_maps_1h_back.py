"""Surface flow maps for medaka 2508 with a 1h backward window.

This is a variant of ``mk2508_surface_flow_maps.py`` that uses a *one-sided*
1h backward temporal window for each key frame:

    window = [frame - WIN_BACK_FRAMES, frame]   # inclusive
            = [t - 60, t]                       # = 1h wide, ending at t

instead of the original ``± HALF_WIN = ± 8`` symmetric window (~16 min
total). The wider 1h window gives a cleaner bulk-motion estimate at each
developmental milestone:

  Frame 116 → window 56..116   = "4 h to 3 h before first contraction"
  Frame 160 → window 100..160  = "1 h before to start of internalisation"
  Frame 236 → window 176..236  = "3 h to 2 h before first contraction"
  Frame 356 → window 296..356  = "2 h to 1 h before first contraction"
  Frame 476 → window 416..476  = "1 h before to first contraction"

Only page 1 (vector flow maps) is produced. There is **no Panel B** in
this output; the existing 2-page PDF from ``mk2508_surface_flow_maps.py``
is preserved untouched.

Outputs (in this directory):
  mk2508_surface_flow_maps_1h_back.pdf   -- single page, vector maps only
  mk2508_surface_flow_maps_1h_back.png   -- same
  mk2508_surface_flow_maps_1h_back_summary.csv
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# Reuse the entire pipeline from the existing module.
import mk2508_surface_flow_maps as m
from mk2508_surface_flow_maps import (
    FRAME_LABELS, KEY_FRAMES,
    VOXEL_XY_UM_PX, VOXEL_Z_UM_PX, FRAME_INTERVAL_SEC, FRAMES_PER_MIN,
    MIN_TRACK_LEN, SPEED_OUTLIER_STD, MAX_INSTANT_SPEED_UM_PER_FRAME,
    FLOW_BIN_UM, FLOW_MIN_N, SURFACE_UM, Z_MID, Z_BAND_UM, EXTENT_PAD_UM,
    ARROW_SCALE, ARROW_LW, ARROW_HEAD_FRAC,
    apply_voxel_calibration, select_tracks, compute_per_step_velocity,
    select_surface, select_slice, bin_xy, classify_bins, draw_panel,
    load_data,
)


# --- 1h-backward window -----------------------------------------------------
WIN_BACK_FRAMES = 60  # = 1h at 30 s/frame
FRAME_LABELS_1H = {
    116: "4 h to 3 h before contraction",
    160: "1 h before to start of internalisation",
    236: "3 h to 2 h before contraction",
    356: "2 h to 1 h before contraction",
    476: "1 h before to first contraction",
}


def select_window_back(
    spots_v: pd.DataFrame,
    frame: int,
    valid_tracks: pd.Series,
    win_back: int = WIN_BACK_FRAMES,
) -> pd.DataFrame:
    """One-sided backward window: [frame - win_back, frame]."""
    lo = frame - win_back
    hi = frame
    return spots_v[
        spots_v["MID_FRAME"].between(lo, hi)
        & spots_v["TRACK_ID"].isin(valid_tracks)
    ]


# Monkey-patch so any downstream helper that calls ``select_window`` from
# the module gets our 1h-back semantics.
m.select_window = select_window_back


HERE = Path(__file__).resolve().parent
IO_DIR = Path("results/medaka_mk2508")
SUFFIX = "1h_back"
STEM = f"mk2508_surface_flow_maps_{SUFFIX}"
OUT_PDF = IO_DIR / f"{STEM}.pdf"
OUT_PNG = IO_DIR / f"{STEM}.png"
OUT_CSV = IO_DIR / f"{STEM}_summary.csv"


def _build_page1(*, spots_v, valid_tracks, sphere, frames,
                 surface_um, z_mid, z_band, bin_um, min_n,
                 arrow_scale, arrow_lw, arrow_head, win_back,
                 n_dropped_speed):
    """Single-page flow-map strip (no Panel B)."""
    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 6.5))
    if n == 1:
        axes = [axes]
    for ax, frame in zip(axes, frames):
        window = select_window_back(spots_v, frame, valid_tracks, win_back)
        surface = select_surface(window, surface_um)
        slice_ = select_slice(surface, z_mid, z_band)
        fl = bin_xy(slice_, bin_um, min_n)
        fl = classify_bins(fl)
        bg = slice_[["MID_X", "MID_Y"]].to_numpy()
        draw_panel(ax, fl, sphere,
                   frame=frame, surface_um=surface_um,
                   z_mid=z_mid, z_band=z_band,
                   bin_um=bin_um, arrow_scale=arrow_scale,
                   arrow_lw=arrow_lw, arrow_head=arrow_head,
                   bg_scatter=bg)
        # Frame label override (1h-window phrasing)
        ax.set_title(
            f"Frame {frame}\n{FRAME_LABELS_1H.get(frame, '')}",
            fontsize=11, fontweight="bold",
        )

    handles = [
        Line2D([0], [0], color="#2166AC", lw=2, label="AP-ward (vy < 0, toward animal pole)"),
        Line2D([0], [0], color="#B2182B", lw=2, label="VP-ward (vy > 0, toward vegetal pole)"),
        Line2D([0], [0], color="black", lw=1, label=f"Sphere (R = {sphere['radius']:.0f} µm)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), frameon=False)

    speed_clause = (
        f"tracks ≥ {MIN_TRACK_LEN} spots; per-step instant speed ≤ "
        f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.1f} µm/frame "
        f"({MAX_INSTANT_SPEED_UM_PER_FRAME * FRAMES_PER_MIN:.0f} µm/min) | "
        f"{n_dropped_speed:,} comet steps dropped"
    )
    fig.suptitle(
        "Medaka 2508 — Surface flow maps at key developmental timepoints "
        "(1h backward window)\n"
        f"Side view (XY, AP at top) | surface band depth ≤ {surface_um:.0f} µm | "
        f"voxel {VOXEL_XY_UM_PX} µm/px | window = [{WIN_BACK_FRAMES} frames, 0] "
        f"({WIN_BACK_FRAMES / FRAMES_PER_MIN:.0f} min) | {speed_clause}",
        fontsize=10.5, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    return fig


def main() -> None:
    io_dir = IO_DIR

    print(f"\n{'='*70}\n  MEDAKA 2508 SURFACE FLOW MAPS — 1h backward window\n{'='*70}\n")

    print("Loading data...")
    spots_px, tracks, sphere_px = load_data(io_dir)
    print(f"  oriented_spots : {len(spots_px):,} rows (raw, in pixels)")
    print(f"  filtered_tracks: {len(tracks):,} rows")
    print(f"  sphere center  : ({sphere_px['center_x']:.1f}, {sphere_px['center_y']:.1f}, "
          f"{sphere_px['center_z']:.1f}) px  R = {sphere_px['radius']:.1f} px")

    print(f"\nVoxel calibration: {VOXEL_XY_UM_PX} µm/px (XY), {VOXEL_Z_UM_PX} µm/px (Z), "
          f"{FRAME_INTERVAL_SEC}s/frame -> {FRAMES_PER_MIN} frames/min")
    spots, sphere = apply_voxel_calibration(spots_px, sphere_px)
    print(f"  sphere R in µm: {sphere['radius']:.1f} µm")

    keep, diag = select_tracks(tracks)
    print("\nTrack filtering (comet / meteor removal):")
    print(f"  total tracks        : {diag['n_total']:,}")
    print(f"  with n_spots >= {MIN_TRACK_LEN:<3} : {diag['n_long']:,}")
    print(f"  tracks removed as comet: {diag['n_comet_removed']:,}")
    print(f"  tracks kept          : {diag['n_kept']:,}  "
          f"({100 * diag['n_kept'] / diag['n_total']:.1f}%)")

    print("\nComputing per-step velocity...")
    spots_v_all = compute_per_step_velocity(spots)
    n_dropped_speed = spots_v_all.attrs.get("n_dropped_speed", 0)
    spots_v = spots_v_all[spots_v_all["TRACK_ID"].isin(keep)]
    print(f"  per-step rows (after kept-track filter)     : {len(spots_v):,}")
    print(f"  per-step rows dropped (instant > "
          f"{MAX_INSTANT_SPEED_UM_PER_FRAME} µm/frame): {n_dropped_speed:,}")

    print(f"\nGenerating flow maps for frames {KEY_FRAMES} "
          f"with 1h backward window [{WIN_BACK_FRAMES} frames, 0]...")
    print(f"  window coverage:")
    for f in KEY_FRAMES:
        lo = f - WIN_BACK_FRAMES
        hi = f
        print(f"    frame {f:>3}: [{lo:>3}, {hi:>3}]  -> {FRAME_LABELS_1H[f]}")

    # ---- Build page 1 figure --------------------------------------------
    fig = _build_page1(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=KEY_FRAMES,
        surface_um=SURFACE_UM, z_mid=Z_MID, z_band=Z_BAND_UM,
        bin_um=FLOW_BIN_UM, min_n=FLOW_MIN_N,
        arrow_scale=ARROW_SCALE, arrow_lw=ARROW_LW,
        arrow_head=ARROW_HEAD_FRAC, win_back=WIN_BACK_FRAMES,
        n_dropped_speed=n_dropped_speed,
    )

    # ---- Save: PDF (single page) + PNG ----------------------------------
    with PdfPages(OUT_PDF) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    # Also a standalone PNG (same content)
    fig2 = _build_page1(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=KEY_FRAMES,
        surface_um=SURFACE_UM, z_mid=Z_MID, z_band=Z_BAND_UM,
        bin_um=FLOW_BIN_UM, min_n=FLOW_MIN_N,
        arrow_scale=ARROW_SCALE, arrow_lw=ARROW_LW,
        arrow_head=ARROW_HEAD_FRAC, win_back=WIN_BACK_FRAMES,
        n_dropped_speed=n_dropped_speed,
    )
    fig2.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    plt.close(fig2)

    # ---- Per-frame summary (binned stats from the 1h-back window) -------
    rows = []
    for f in KEY_FRAMES:
        window = select_window_back(spots_v, f, keep, WIN_BACK_FRAMES)
        surface = select_surface(window, SURFACE_UM)
        slice_ = select_slice(surface, Z_MID, Z_BAND_UM)
        fl = bin_xy(slice_, FLOW_BIN_UM, FLOW_MIN_N)
        fl = classify_bins(fl)
        if len(fl) > 0:
            rows.append({
                "frame": f,
                "label_1h_window": FRAME_LABELS_1H[f],
                "label_original": FRAME_LABELS[f],
                "win_lo": f - WIN_BACK_FRAMES,
                "win_hi": f,
                "n_vectors": int(len(slice_)),
                "n_bins": int(len(fl)),
                "n_ap_ward_bins": int((fl["move_type"] == "AP-ward").sum()),
                "frac_ap_ward_bins": float((fl["move_type"] == "AP-ward").mean()),
                "mean_bin_speed_um_per_frame": float(fl["speed"].mean()),
                "mean_bin_speed_um_per_min": float(fl["speed"].mean() * FRAMES_PER_MIN),
                "mean_bin_vx_um_per_frame": float(fl["vx"].mean()),
                "mean_bin_vy_um_per_frame": float(fl["vy"].mean()),
            })
        else:
            rows.append({
                "frame": f,
                "label_1h_window": FRAME_LABELS_1H[f],
                "label_original": FRAME_LABELS[f],
                "win_lo": f - WIN_BACK_FRAMES,
                "win_hi": f,
                "n_vectors": 0, "n_bins": 0, "n_ap_ward_bins": 0,
                "frac_ap_ward_bins": np.nan,
                "mean_bin_speed_um_per_frame": np.nan,
                "mean_bin_speed_um_per_min": np.nan,
                "mean_bin_vx_um_per_frame": np.nan,
                "mean_bin_vy_um_per_frame": np.nan,
            })
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUT_CSV, index=False)

    print("\nPer-frame summary (1h backward window):")
    print(summary_df.to_string(index=False))
    print(f"\nSaved PDF to {OUT_PDF}")
    print(f"Saved PNG to {OUT_PNG}")
    print(f"Saved summary to {OUT_CSV}")


if __name__ == "__main__":
    main()
