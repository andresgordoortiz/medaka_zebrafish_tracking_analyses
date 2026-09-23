"""Surface flow maps for medaka 2805 minimal with a 1h backward window.

This is the **2805 minimal** counterpart of
``results/medaka_mk2508/mk2508_surface_flow_maps_1h_back.py``.

Key adaptations vs. the 2508 script:

  - The 2805 dataset ships positions **already in µm** (not pixels), so
    ``apply_voxel_calibration`` is an identity transform (voxel = 1.0 µm/px).
  - There is no separate ``filtered_tracks.csv``; we derive track-level
    ``NUMBER_SPOTS`` and ``TRACK_MEAN_SPEED`` directly from
    ``oriented_spots.csv`` via groupby.
  - The dataset ends at frame 475 (vs. 697 for 2508), so the last key frame
    is 475 (not 476). The other anchors (116, 160, 236, 356) are inherited
    from 2508 — the user can override via ``--frames``.
  - Median track length is only 5 spots (vs. 10 for 2508), so we drop
    ``MIN_TRACK_LEN`` to 5 by default. The per-step speed cap of
    ``MAX_INSTANT_SPEED_UM_PER_FRAME = 3.0 µm/frame`` is kept.

Outputs (in this directory):
  mk2805_surface_flow_maps_1h_back.pdf       -- single page, vector maps only
  mk2805_surface_flow_maps_1h_back.png       -- same
  mk2805_surface_flow_maps_1h_back_summary.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Locate the 2508 pipeline source we will reuse for the heavy lifting.
# The 2508 script lives in the same directory as this 2805 script.
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
M2508_DIR = HERE
sys.path.insert(0, str(M2508_DIR))

import mk2508_surface_flow_maps as m  # noqa: E402
from mk2508_surface_flow_maps import (  # noqa: E402
    KEY_FRAMES, FRAME_LABELS,
    VOXEL_XY_UM_PX, VOXEL_Z_UM_PX, FRAME_INTERVAL_SEC, FRAMES_PER_MIN,
    MIN_TRACK_LEN, SPEED_OUTLIER_STD, MAX_INSTANT_SPEED_UM_PER_FRAME,
    FLOW_BIN_UM, FLOW_MIN_N, SURFACE_UM, Z_MID, Z_BAND_UM, EXTENT_PAD_UM,
    ARROW_SCALE, ARROW_LW, ARROW_HEAD_FRAC,
    compute_per_step_velocity,
    select_surface, select_slice, bin_xy, classify_bins, draw_panel,
)

# ---------------------------------------------------------------------------
# 2805-specific overrides
# ---------------------------------------------------------------------------
# Positions are already in µm; the calibration step is an identity.
VOXEL_XY_UM_PX_2805 = 1.0
VOXEL_Z_UM_PX_2805 = 1.0

# Last key frame is 475 (vs. 476 in 2508, which has frames 0..697).
KEY_FRAMES_2805 = [116, 160, 236, 356, 475]

# Track filter: median track length is only 5 in this dataset, so relax to 5.
MIN_TRACK_LEN_2805 = 5

# 1h backward window at 30 s/frame.
WIN_BACK_FRAMES = 60  # = 1h at 30 s/frame

# Frame labels for the 1h-backward window phrasing. Kept generic ("h before
# contraction") since the 2805 developmental timeline has not been anchored
# yet relative to 2508 — user can refine later.
FRAME_LABELS_2805 = {
    116: "3 h before first contraction",
    160: "Start of internalisation",
    236: "2 h before first contraction",
    356: "Mid-epiboly / pre-contraction",
    475: "First contraction (end)",
}
FRAME_LABELS_1H = {
    116: "4 h to 3 h before contraction",
    160: "1 h before to start of internalisation",
    236: "3 h to 2 h before contraction",
    356: "2 h to 1 h before contraction",
    475: "1 h before to first contraction",
}


# ---------------------------------------------------------------------------
# I/O overrides for the 2805 minimal dataset
# ---------------------------------------------------------------------------
def load_data_2805(io_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return (oriented_spots in µm, derived track summary in µm, sphere).

    The 2805 minimal dataset has:
      - oriented_spots.csv   — already in µm
      - oriented_tracks.csv  — same data as oriented_spots.csv
      - sphere_params.csv    — center and R already in µm

    We build the per-track summary (NUMBER_SPOTS, TRACK_MEAN_SPEED) by
    groupby on the spots, mimicking the columns that
    ``select_tracks`` expects from ``filtered_tracks.csv``.
    """
    spots = pd.read_csv(io_dir / "oriented_spots.csv", low_memory=False)
    # Per-track stats from the per-step velocity table (vx, vy in µm/frame).
    # We compute this lazily; for now just NUMBER_SPOTS, TRACK_MEAN_SPEED
    # in displacement-per-frame units.
    g = spots.groupby("TRACK_ID", sort=False)
    track_summary = pd.DataFrame({
        "TRACK_ID":      [k for k, _ in g],
        "NUMBER_SPOTS":  g.size().to_numpy(),
    })
    # Median instantaneous speed in µm/frame per track = mean of sqrt(dx^2 + dy^2)/dt
    # over each track's steps.
    sorted_sp = spots.sort_values(["TRACK_ID", "FRAME"]).copy()
    sorted_sp["dx"] = sorted_sp.groupby("TRACK_ID")["POSITION_X"].diff()
    sorted_sp["dy"] = sorted_sp.groupby("TRACK_ID")["POSITION_Y"].diff()
    sorted_sp["dt"] = sorted_sp.groupby("TRACK_ID")["FRAME"].diff()
    sorted_sp["step_speed"] = (
        np.sqrt(sorted_sp["dx"] ** 2 + sorted_sp["dy"] ** 2) / sorted_sp["dt"]
    )
    track_summary = track_summary.merge(
        sorted_sp.groupby("TRACK_ID")["step_speed"].mean()
        .rename("TRACK_MEAN_SPEED").reset_index(),
        on="TRACK_ID",
    )

    sph = pd.read_csv(io_dir / "sphere_params.csv")
    sphere = {row["parameter"]: row["value"] for _, row in sph.iterrows()}
    for k in ("radius", "center_x", "center_y", "center_z",
              "rmse", "cap_height", "cap_base_radius"):
        if k in sphere:
            try:
                sphere[k] = float(sphere[k])
            except (TypeError, ValueError):
                pass
    return spots, track_summary, sphere


def apply_voxel_calibration_2805(spots: pd.DataFrame, sphere: dict
                                 ) -> tuple[pd.DataFrame, dict]:
    """Identity transform for positions already in µm."""
    return spots.copy(), sphere.copy()


def select_tracks_2805(tracks: pd.DataFrame,
                       min_len: int = MIN_TRACK_LEN_2805,
                       k_std: float = SPEED_OUTLIER_STD
                       ) -> tuple[pd.Series, dict]:
    """Return (mask of valid track IDs, diagnostics dict) — µ m/frame units."""
    n0 = tracks["TRACK_ID"].nunique()
    long_enough = tracks["NUMBER_SPOTS"] >= min_len
    speeds = tracks["TRACK_MEAN_SPEED"].astype(float)
    median = float(speeds.median())
    std = float(speeds.std())
    thr = median + k_std * std
    not_comet = speeds <= thr
    keep = tracks[long_enough & not_comet]["TRACK_ID"]
    diag = {
        "n_total": n0,
        "n_long": int(long_enough.sum()),
        "median_speed": median,
        "std_speed": std,
        "speed_thr": thr,
        "n_comet_removed": int((long_enough & ~not_comet).sum()),
        "n_kept": int(len(keep)),
    }
    return keep, diag


def select_window_back_2805(
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


# Monkey-patch the 2508 module so any helper that calls these by name gets
# 2805 behaviour.  This keeps the rest of the pipeline untouched.
m.load_data = load_data_2805
m.apply_voxel_calibration = apply_voxel_calibration_2805
m.select_tracks = select_tracks_2805
m.select_window = select_window_back_2805


# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
SUFFIX = "1h_back"
STEM = f"mk2805_surface_flow_maps_{SUFFIX}"
OUT_PDF = Path("results/medaka_mk2805_minimal") / f"{STEM}.pdf"
OUT_PNG = Path("results/medaka_mk2805_minimal") / f"{STEM}.png"
OUT_CSV = Path("results/medaka_mk2805_minimal") / f"{STEM}_summary.csv"
IO_DIR  = Path("results/medaka_mk2805_minimal")


# ---------------------------------------------------------------------------
# Figure builder
# ---------------------------------------------------------------------------
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
        window = select_window_back_2805(spots_v, frame, valid_tracks, win_back)
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
        ax.set_title(
            f"Frame {frame}\n{FRAME_LABELS_1H.get(frame, '')}",
            fontsize=11, fontweight="bold",
        )

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
        f"tracks ≥ {MIN_TRACK_LEN_2805} spots; per-step instant speed ≤ "
        f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.1f} µm/frame "
        f"({MAX_INSTANT_SPEED_UM_PER_FRAME * FRAMES_PER_MIN:.0f} µm/min) | "
        f"{n_dropped_speed:,} comet steps dropped"
    )
    fig.suptitle(
        "Medaka 2805 — Surface flow maps at key developmental timepoints "
        "(1h backward window)\n"
        f"Side view (XY, AP at top) | surface band depth ≤ {surface_um:.0f} µm | "
        f"voxel 1.0 µm/px (positions already in µm) | "
        f"window = [{WIN_BACK_FRAMES} frames, 0] "
        f"({WIN_BACK_FRAMES / FRAMES_PER_MIN:.0f} min) | {speed_clause}",
        fontsize=10.5, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    io_dir = IO_DIR

    print(f"\n{'='*70}\n  MEDAKA 2805 SURFACE FLOW MAPS — 1h backward window\n{'='*70}\n")

    print("Loading data...")
    spots_um, tracks_um, sphere_um = load_data_2805(io_dir)
    print(f"  oriented_spots : {len(spots_um):,} rows (already in µm)")
    print(f"  derived tracks : {len(tracks_um):,} tracks")
    print(f"  sphere center  : ({sphere_um['center_x']:.1f}, "
          f"{sphere_um['center_y']:.1f}, {sphere_um['center_z']:.1f}) µm  "
          f"R = {sphere_um['radius']:.1f} µm")

    # Skip voxel calibration (positions already in µm).
    spots, sphere = apply_voxel_calibration_2805(spots_um, sphere_um)

    keep, diag = select_tracks_2805(tracks_um)
    print("\nTrack filtering (comet / meteor removal):")
    print(f"  total tracks         : {diag['n_total']:,}")
    print(f"  with n_spots >= {MIN_TRACK_LEN_2805:<3} : {diag['n_long']:,}")
    print(f"  TRACK_MEAN_SPEED median = {diag['median_speed']:.4f}  "
          f"std = {diag['std_speed']:.4f} (µm/frame)")
    print(f"  threshold = median + {SPEED_OUTLIER_STD}·σ = {diag['speed_thr']:.4f} µm/frame "
          f"({diag['speed_thr'] * FRAMES_PER_MIN:.4f} µm/min)")
    print(f"  tracks removed as comet: {diag['n_comet_removed']:,}")
    print(f"  tracks kept           : {diag['n_kept']:,}  "
          f"({100 * diag['n_kept'] / diag['n_total']:.1f}%)")

    print("\nComputing per-step velocity...")
    spots_v_all = compute_per_step_velocity(spots)
    n_dropped_speed = spots_v_all.attrs.get("n_dropped_speed", 0)
    spots_v = spots_v_all[spots_v_all["TRACK_ID"].isin(keep)]
    print(f"  per-step rows (after kept-track filter)     : {len(spots_v):,}")
    print(f"  per-step rows dropped (instant > "
          f"{MAX_INSTANT_SPEED_UM_PER_FRAME} µm/frame): {n_dropped_speed:,}")
    print(f"  instantaneous speed stats (µm / frame):")
    print(f"    median = {spots_v['speed'].median():.3f}   "
          f"mean = {spots_v['speed'].mean():.3f}   "
          f"std = {spots_v['speed'].std():.3f}")
    print(f"    99th pct = {spots_v['speed'].quantile(0.99):.3f}   "
          f"max = {spots_v['speed'].max():.3f}")

    print(f"\nGenerating flow maps for frames {KEY_FRAMES_2805} "
          f"with 1h backward window [{WIN_BACK_FRAMES} frames, 0]...")
    print(f"  window coverage:")
    for f in KEY_FRAMES_2805:
        lo = f - WIN_BACK_FRAMES
        hi = f
        print(f"    frame {f:>3}: [{lo:>3}, {hi:>3}]  -> {FRAME_LABELS_1H[f]}")

    # ---- Build page 1 figure ---------------------------------------------
    fig = _build_page1(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=KEY_FRAMES_2805,
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
    fig2 = _build_page1(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=KEY_FRAMES_2805,
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
    for f in KEY_FRAMES_2805:
        window = select_window_back_2805(spots_v, f, keep, WIN_BACK_FRAMES)
        surface = select_surface(window, SURFACE_UM)
        slice_ = select_slice(surface, Z_MID, Z_BAND_UM)
        fl = bin_xy(slice_, FLOW_BIN_UM, FLOW_MIN_N)
        fl = classify_bins(fl)
        if len(fl) > 0:
            rows.append({
                "frame": f,
                "label_1h_window": FRAME_LABELS_1H[f],
                "label_original": FRAME_LABELS_2805[f],
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
                "label_original": FRAME_LABELS_2805[f],
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