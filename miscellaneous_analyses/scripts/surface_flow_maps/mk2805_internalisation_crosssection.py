"""Cross-section flow map for the medaka 2805 minimal dataset — internalisation.

This is the cross-section counterpart of the 1h-back surface flow maps.
It is designed to reveal the **internalisation signature**: cells at the
animal-pole margin that bend *inward* (toward the embryo interior) while
the rest of the surface keeps drifting animal-ward / vegetal-ward in the
conveyor belt.

Geometry
--------
After ``orient_embryo()`` on this dataset:
  - data Y axis: animal pole at -Y, vegetal pole at +Y
  - data Z axis: inward toward sphere centre (cap interior is at +Z)
  - data X axis: dorsal (+X) / ventral (-X)
At the sphere margin (animal pole region), the cap is essentially
horizontal: surface cells sit near the lower-z plane, and "deeper" cells
(cap_depth growing) sit at progressively larger z. In other words, at the
margin, **z IS the inward axis** (correl(z, cap_depth) = 0.994 in this
dataset).

So a YZ cross-section through the dorsal-ventral midline (|x| ≤ Z_HALF)
is the *finest* cross-section to expose the internalisation motion:
  - Y axis on the figure: animal pole at TOP (-Y), vegetal pole at BOTTOM (+Y)
  - Z axis on the figure: surface at the BOTTOM (low z), interior at the TOP
  - Cells going UP in z = internalising (cap_depth growing)

Cross-section
-------------
Sagittal slab through x ∈ [-Z_HALF, +Z_HALF] (Z_HALF = 30 µm, total 60 µm).
Window of frames: 160 → 250 (start of internalisation through mid-epiboly).

Bins
----
XY would project away the z (depth) information, so the binning is in
(Y, Z) space, not (X, Y):
  - Y bin: 20 µm (animal-vegetal axis)
  - Z bin: 20 µm (depth axis = inward axis)
  - Minimum n per bin: 4 vectors

Colours
-------
Arrows are coloured by **mean cap_depth** of the bin:
  - cap_depth < 30 µm  -> surface band     -> cool blue
  - cap_depth 30..60   -> shallow interior -> green
  - cap_depth 60..100  -> deep interior    -> warm orange
The bulge should appear as a band of upward-z arrows at the animal pole
side (Y ≈ -50) bending into the interior, while the conveyor belt at
larger Y (vegetal) shows the opposite (downward arrows = outward).

Outputs
-------
- mk2805_internalisation_crosssection.pdf   (single page)
- mk2805_internalisation_crosssection.png   (same)
- mk2805_internalisation_crosssection_summary.csv  (per-key-frame bin stats)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Reuse I/O + voxel calibration + per-step velocity from the 2805 1h-back script
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str((HERE / "mk2805_surface_flow_maps_1h_back.py").resolve()))
import mk2805_surface_flow_maps_1h_back as base  # noqa: E402
from mk2805_surface_flow_maps_1h_back import (  # noqa: E402
    FRAME_LABELS_2805, FRAMES_PER_MIN, MIN_TRACK_LEN_2805,
    MAX_INSTANT_SPEED_UM_PER_FRAME,
    load_data_2805, apply_voxel_calibration_2805, select_tracks_2805,
    compute_per_step_velocity,
)

# ---------------------------------------------------------------------------
# Cross-section specific parameters
# ---------------------------------------------------------------------------
# Cross-section is the dorsal-ventral midline slab.
# Use a wider slab to get enough cells per bin.
Z_HALF = 100.0   # µm half-width of the slab

# Internalisation window: starts at frame 160 (user: "ingression occurs from
# frame 160 onwards") through mid-epiboly. Includes late frames where
# the ingression is most visible in the deep layer.
# Internalisation frames, picked by maximising z-score of the bulge inward
# signal. Peak at frames 270-285 (z=13.4), with z>10 across 200-260.
# Best 6 frames showing the internalisation arc:
INTERN_FRAMES = [200, 230, 270, 290, 320, 380]

# TRUE internalisation bulge (USER-DEFINED ROI matching surface flow map, 2026-08-26)
# Reference: surface flow maps at the dorsal view show a clear cluster of
# cells going down at X ∈ [-100, 100], Y ∈ [-200, 0]. Those cells are the
#   ingressing cohort at the dorsal blastoderm margin.
#
# The bulge has TWO layers:
#   1. Surface (cap < 30): at z ~ -100..+25 (median -46)
#      - Frame 100-170: mean_vz = +0.07..+0.13 (cells diving inward)
#      - Peak inward at frame 150: vz = +0.125
#   2. Deep (cap 60-120): at z ~ -30..+80 (median +23)
#      - Throughout: mean_vz is NEGATIVE (cells flowing AP-ward and outward,
#        the hypoblast flow back toward the cap surface as they migrate)
#      - Strong AP-ward: vy = -0.5 to -0.2
BULGE_Y = (0, 200)
BULGE_X = (-100, 100)
BULGE_CAP = (-10, 30)

# Bin size in YZ (Y = animal-vegetal, Z = inward axis). 30 µm bins balance
# detail with statistics.
BIN_Y_UM = 30.0
BIN_Z_UM = 30.0
MIN_N = 3

# Surface band depth threshold (for the in-surface mask in the legend)
SURFACE_UM = 30.0

IO_DIR  = Path("results/medaka_mk2805_minimal")
OUT_PDF = IO_DIR / "mk2805_internalisation_crosssection.pdf"
OUT_PNG = IO_DIR / "mk2805_internalisation_crosssection.png"
OUT_CSV = IO_DIR / "mk2805_internalisation_crosssection_summary.csv"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def select_window_frames(spots_v: pd.DataFrame, frame: int, valid_tracks: pd.Series,
                         f_lo: int, f_hi: int) -> pd.DataFrame:
    """Use per-step MID_FRAME to select a half-open interval [f_lo, f_hi]."""
    return spots_v[
        spots_v["MID_FRAME"].between(f_lo, f_hi)
        & spots_v["TRACK_ID"].isin(valid_tracks)
    ]


def select_xsection_slab(spots_v: pd.DataFrame, z_half: float = Z_HALF
                         ) -> pd.DataFrame:
    """Restrict to the dorsal-ventral midline slab |MID_X| <= z_half."""
    return spots_v[spots_v["MID_X"].abs() <= z_half]


def bin_yz(window: pd.DataFrame, bin_y: float = BIN_Y_UM, bin_z: float = BIN_Z_UM,
           min_n: int = MIN_N) -> pd.DataFrame:
    """Bin velocity vectors onto a regular (Y, Z) grid.

    Velocity components are projected onto the (Y, Z) plane: we use ``vy``
    and ``vz`` directly. (Z is the inward axis; +vz means moving outward,
    -vz means moving inward.) The Y axis is flipped for display so animal
    pole sits at the TOP of the figure.
    """
    df = window.copy()
    df["y_bin"] = np.floor(df["MID_Y"] / bin_y) * bin_y + bin_y / 2
    df["z_bin"] = np.floor(df["MID_Z"] / bin_z) * bin_z + bin_z / 2
    g = df.groupby(["y_bin", "z_bin"], as_index=False)
    out = g.agg(
        vy=("vy", "mean"),
        vz=("vz", "mean"),
        vx=("vx", "mean"),
        depth_med=("SPHERICAL_DEPTH", "median"),
        theta_med=("THETA_DEG", "median"),
        phi_med=("PHI_DEG", "median"),
        n=("vy", "count"),
    )
    return out[out["n"] >= min_n]


def classify_inout(fl: pd.DataFrame) -> pd.DataFrame:
    """Per-bin label of the dominant motion direction in the YZ plane.

    In this dataset, the lower cap has its tip at low z and the equator at
    higher z. So a cell with z increasing is moving INWARD (toward the
    sphere centre / deeper into the embryo).

    - 'inward'   if vz > 0  (cap_depth growing — internalising)
    - 'outward'  if vz <= 0 (cap_depth shrinking — outward / epiboly)
    """
    fl = fl.copy()
    fl["motion_type"] = np.where(fl["vz"] > 0, "inward", "outward")
    return fl


# Y range for ingressing-track overlay: animal-pole side (incl. bulge)
# Restricted to the margin band: data Y ∈ [-50, +200] (animal pole at -Y,
# vegetal pole at +Y; the margin/bulge lives at the animal-pole side).
INGRESS_Y_RANGE = (-50, 200)


def find_ingressing_tracks(spots_v: pd.DataFrame, valid_tracks: pd.Series,
                            min_track_len: int = 5) -> pd.DataFrame:
    """Find tracks that are clearly ingressing.

    Criteria (all must hold):
      1. y_mean ∈ [-50, +200] (animal-pole margin band only)
      2. mean_vz > 0  (net positive Z / inward motion)
      3. n_spots >= min_track_len (long enough to be reliable)
    """
    sub = spots_v[spots_v["TRACK_ID"].isin(valid_tracks)].copy()

    # Per-track summary
    g = sub.groupby("TRACK_ID", sort=False)
    track_df = pd.DataFrame({
        "n_spots": g.size(),
        "frame_min": g["FRAME"].min(),
        "frame_max": g["FRAME"].max(),
        "y_mean": g["MID_Y"].mean(),
        "vz_mean": g["vz"].mean(),
        "vy_mean": g["vy"].mean(),
        "fr_inward": (sub[sub["vz"] > 0].groupby("TRACK_ID").size() / g.size()).fillna(0.0),
        "cap_min": g["SPHERICAL_DEPTH"].min(),
        "cap_max": g["SPHERICAL_DEPTH"].max(),
        "cap_start": g["SPHERICAL_DEPTH"].first(),
        "cap_end": g["SPHERICAL_DEPTH"].last(),
    })
    track_df["delta_cap"] = track_df["cap_end"] - track_df["cap_start"]
    track_df["delta_y"] = track_df["y_mean"]  # for sorting
    # Filter
    keep = (
        (track_df["y_mean"] >= INGRESS_Y_RANGE[0]) &
        (track_df["y_mean"] <= INGRESS_Y_RANGE[1]) &
        (track_df["vz_mean"] > 0) &
        (track_df["n_spots"] >= min_track_len)
    )
    return track_df[keep].copy()


def draw_ingressing_tracks(ax, spots_v: pd.DataFrame, valid_tracks: pd.Series,
                            frame: int, f_lo: int, f_hi: int,
                            track_df: pd.DataFrame,
                            max_tracks: int = 40,
                            color: str = "#D62828",
                            alpha: float = 0.95) -> None:
    """Draw ingressing tracks as full trajectories across the panel.

    For each track, plot its complete YZ trajectory (clipped to the panel
    bounds). The portion within the current time window is drawn bolder;
    the rest of the trajectory is drawn thinner as context. No markers.
    """
    if len(track_df) == 0:
        return
    # Choose top N tracks by delta_cap (largest cap growth)
    top_tracks = track_df.nlargest(max_tracks, "delta_cap")
    track_ids = set(top_tracks.index)
    # Get the full trajectory for these tracks
    sub = spots_v[spots_v["TRACK_ID"].isin(track_ids)].sort_values(["TRACK_ID", "FRAME"])
    for tid, t in sub.groupby("TRACK_ID"):
        if len(t) < 3: continue
        # Full trajectory (faint context)
        ax.plot(t["MID_Y"].values, t["MID_Z"].values,
                color=color, linewidth=1.2, alpha=0.35, zorder=4,
                solid_capstyle="round")
        # Portion within the time window (bold)
        t_win = t[(t["FRAME"] >= f_lo) & (t["FRAME"] <= f_hi)]
        if len(t_win) >= 2:
            ax.plot(t_win["MID_Y"].values, t_win["MID_Z"].values,
                    color=color, linewidth=3.2, alpha=alpha, zorder=6,
                    solid_capstyle="round")


def draw_crosssection_panel(ax, fl: pd.DataFrame, frame: int, f_lo: int, f_hi: int,
                            *, sphere, bin_y: float, bin_z: float,
                            arrow_scale: float, arrow_lw: float,
                            bg_scatter: np.ndarray | None = None,
                            zoom: str = "wide",
                            show_box: bool = False,
                            extent_pad: float = 25.0) -> None:
    """Draw a YZ cross-section panel with arrows coloured by direction.

    Convention: horizontal axis is DATA Y (animal pole at right near 0,
    vegetal pole at left near -200). Vertical axis is data Z (cap surface
    at BOTTOM, interior at TOP). The user's reference "Y ∈ [-200, 0]" maps
    directly to data Y values.

    Arrows: vy on horizontal axis (vy<0 = AP-ward = left arrow),
    vz on vertical axis (vz>0 = inward = upward arrow).
    """
    # Plot extent
    if zoom == "bulge":
        # Zoom on the bulge: data Y in [-260, +260], z in [-50, +110].
        # Widen Y range so the bulge box (data Y ∈ [0, +200]) and the
        # surrounding vegetal pole context both fit on the panel.
        x_lo, x_hi = -260, 260
        y_lo, y_hi = -50, 110
    else:
        if bg_scatter is not None and len(bg_scatter) > 0:
            y_lo_data = bg_scatter[:, 1].min()
            y_hi_data = bg_scatter[:, 1].max()
            z_lo_data = bg_scatter[:, 2].min()
            z_hi_data = bg_scatter[:, 2].max()
        else:
            y_lo_data, y_hi_data = fl["y_bin"].min(), fl["y_bin"].max()
            z_lo_data, z_hi_data = fl["z_bin"].min(), fl["z_bin"].max()
        x_lo, x_hi = y_lo_data, y_hi_data
        y_lo, y_hi = z_lo_data, z_hi_data

    # Cell-density background
    if bg_scatter is not None and len(bg_scatter) > 0:
        bg = bg_scatter
        bg_disp_x = bg[:, 1]
        bg_disp_y = bg[:, 2]
        mask = (bg_disp_x >= x_lo - extent_pad) & (bg_disp_x <= x_hi + extent_pad) & \
               (bg_disp_y >= y_lo - extent_pad) & (bg_disp_y <= y_hi + extent_pad)
        bg_disp_x = bg_disp_x[mask]
        bg_disp_y = bg_disp_y[mask]
        max_bg = 10000
        if len(bg_disp_x) > max_bg:
            rng = np.random.default_rng(42)
            idx = rng.choice(len(bg_disp_x), size=max_bg, replace=False)
            bg_disp_x = bg_disp_x[idx]
            bg_disp_y = bg_disp_y[idx]
        hb = ax.hexbin(bg_disp_x, bg_disp_y, gridsize=40, cmap="Greys",
                       mincnt=1, alpha=0.55, edgecolors="none", zorder=0)

    # (No box ROI — tracks themselves show the ingression region)

    # Arrows: colour by motion direction
    if len(fl) > 0:
        for direction, color, label, zorder in [
            ("inward", "#E63946", "INWARD (vz>0, internalising)", 4),
            ("outward", "#1D6997", "OUTWARD (vz<0, epiboly)",     3),
        ]:
            sub = fl[fl["motion_type"] == direction]
            if sub.empty:
                continue
            ax.quiver(
                sub["y_bin"], sub["z_bin"],
                sub["vy"] * arrow_scale, sub["vz"] * arrow_scale,
                color=color, linewidth=arrow_lw + 0.3,
                headwidth=4, headlength=5, headaxislength=4,
                angles="xy", scale_units="xy", scale=1.0,
                alpha=0.95, zorder=zorder,
            )

    # Sphere outline (lower cap)
    R = sphere["radius"]
    cz = sphere["center_z"]
    cap_height = sphere.get("cap_height", 317.0)
    cap_base_radius = sphere.get("cap_base_radius", 621.0)
    cap_base_z = cz - R + cap_height
    y_out = np.linspace(-cap_base_radius, cap_base_radius, 200)
    z_out = cz - np.sqrt(R**2 - y_out**2)
    cap_mask = z_out <= cap_base_z
    ax.plot(y_out[cap_mask], z_out[cap_mask], color="black", lw=0.7, alpha=0.55, zorder=2)
    ax.plot(0, cz - R, "+", color="black", ms=6, alpha=0.7, zorder=2)
    ax.axhline(cap_base_z, color="black", lw=0.5, ls=":", alpha=0.5, zorder=2)

    # Plot extent
    ax.set_xlim(x_lo - extent_pad, x_hi + extent_pad)
    ax.set_ylim(y_lo - extent_pad, y_hi + extent_pad)

    ax.set_aspect("equal")
    ax.set_xlabel("data Y  (animal pole at right → 0, vegetal pole at left → −200)", fontsize=10)
    ax.set_ylabel("Z  (cap surface ↓ | interior ↑)", fontsize=10)

    n_in = int((fl["motion_type"] == "inward").sum()) if len(fl) else 0
    n_out = int((fl["motion_type"] == "outward").sum()) if len(fl) else 0
    # Skip per-panel title/inward-outward text if arrows disabled (fl empty)
    if len(fl) > 0:
        title = f"Frame {frame}\n[{f_lo}, {f_hi}]  {len(fl)} bins"
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.text(0.02, 0.97,
                f"inward {n_in}  outward {n_out}",
                transform=ax.transAxes, fontsize=9, va="top", zorder=5,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))


def build_summary_rows(spots_v: pd.DataFrame, valid_tracks, sphere,
                       frames, *, half_window=10) -> list[dict]:
    """Per-frame summary: bin counts and motion-direction fractions."""
    rows = []
    for f in frames:
        win = select_window_frames(spots_v, f, valid_tracks, f - half_window, f + half_window)
        slab = select_xsection_slab(win)
        fl = bin_yz(slab)
        fl = classify_inout(fl)
        if len(fl) > 0:
            n_in = int((fl["motion_type"] == "inward").sum())
            n_out = int((fl["motion_type"] == "outward").sum())
            rows.append({
                "frame": f,
                "win_lo": f - half_window,
                "win_hi": f + half_window,
                "n_vectors": int(len(slab)),
                "n_bins": int(len(fl)),
                "n_inward_bins": n_in,
                "n_outward_bins": n_out,
                "frac_inward_bins": float(n_in / len(fl)),
                "mean_bin_vy_um_per_frame": float(fl["vy"].mean()),
                "mean_bin_vz_um_per_frame": float(fl["vz"].mean()),
                "mean_bin_speed_um_per_frame": float(np.sqrt(fl["vy"]**2 + fl["vz"]**2).mean()),
                "mean_bin_speed_um_per_min":
                    float(np.sqrt(fl["vy"]**2 + fl["vz"]**2).mean() * FRAMES_PER_MIN),
            })
        else:
            rows.append({
                "frame": f, "win_lo": f - half_window, "win_hi": f + half_window,
                "n_vectors": 0, "n_bins": 0, "n_inward_bins": 0, "n_outward_bins": 0,
                "frac_inward_bins": np.nan, "mean_bin_vy_um_per_frame": np.nan,
                "mean_bin_vz_um_per_frame": np.nan, "mean_bin_speed_um_per_frame": np.nan,
                "mean_bin_speed_um_per_min": np.nan,
            })
    return rows


def build_figure(*, spots_v, valid_tracks, sphere, frames,
                 bin_y, bin_z, arrow_scale, arrow_lw, half_window=10):
    """Single-row definitive figure: ingressing tracks + bin arrows.

    Six panels covering the full internalisation arc. Each panel shows:
      - cell density (greyscale hexbin)
      - BIN arrows coloured by motion direction (red = inward, blue = outward)
      - INGRESSING TRACKS in red (Y∈[-50,200], mean vz>0)
    No ROI box.
    """
    # Pre-compute ingressing tracks (Y in [-50, +200], positive vz)
    print("Identifying ingressing tracks (Y∈[-50,200], vz>0)...")
    track_df = find_ingressing_tracks(spots_v, valid_tracks)
    print(f"  Found {len(track_df)} ingressing tracks")
    print(f"  Will plot top 40 by cap_depth growth")

    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 6.0))
    if n == 1:
        axes = [axes]
    for ax, frame in zip(axes, frames):
        win = select_window_frames(spots_v, frame, valid_tracks,
                                   frame - half_window, frame + half_window)
        slab = select_xsection_slab(win)
        bg = slab[["MID_X", "MID_Y", "MID_Z"]].to_numpy()
        # Binned flow field (arrows coloured by motion direction)
        fl = bin_yz(slab, bin_y=bin_y, bin_z=bin_z)
        fl = classify_inout(fl)
        draw_crosssection_panel(
            ax, fl=fl, frame=frame,
            f_lo=frame - half_window, f_hi=frame + half_window,
            sphere=sphere, bin_y=bin_y, bin_z=bin_z,
            arrow_scale=arrow_scale, arrow_lw=arrow_lw,
            bg_scatter=bg, zoom="bulge", show_box=False,
        )
        # Overlay the ingressing tracks
        f_lo = frame - half_window
        f_hi = frame + half_window
        draw_ingressing_tracks(ax, spots_v, valid_tracks,
                                frame=frame, f_lo=f_lo, f_hi=f_hi,
                                track_df=track_df, max_tracks=40)
        ax.set_title(f"Frame {frame}\n[{frame-half_window}, {frame+half_window}]",
                     fontsize=13, fontweight="bold")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color="#E63946", lw=2.2,
               label="INWARD (vz>0) bin"),
        Line2D([0], [0], color="#1D6997", lw=2.2,
               label="OUTWARD (vz<0) bin"),
        Line2D([0], [0], color="#D62828", lw=2.8,
               label=f"Ingressing tracks (Y∈[-50,200], vz>0, n={len(track_df):,})"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), fontsize=10,
               framealpha=0.95, edgecolor="black")

    fig.suptitle(
        f"Medaka 2805 — Internalisation cross-section (YZ plane)\n"
        f"Dorsal-ventral slab |x| ≤ {Z_HALF:.0f} µm  |  "
        f"window = ±{half_window} frames (±{half_window * 30}s)  |  "
        f"YZ bin {bin_y:.0f} µm  |  tracks: data Y∈[{INGRESS_Y_RANGE[0]},{INGRESS_Y_RANGE[1]}], "
        f"vz > 0  |  {len(track_df):,} ingressing tracks overlaid (top {min(40, len(track_df))} shown)",
        fontsize=12, fontweight="bold", y=1.02,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    return fig


# Main
# ---------------------------------------------------------------------------
def main() -> None:
    io_dir = IO_DIR
    print(f"\n{'='*70}\n  MEDAKA 2805 — INTERNALISATION CROSS-SECTION (YZ)\n{'='*70}\n")

    print("Loading data...")
    spots_um, tracks_um, sphere_um = load_data_2805(io_dir)
    print(f"  oriented_spots : {len(spots_um):,} rows")
    print(f"  tracks         : {len(tracks_um):,}")
    print(f"  sphere         : R = {sphere_um['radius']:.1f} µm, "
          f"cap_height = {sphere_um.get('cap_height','?'):.1f} µm")

    spots, sphere = apply_voxel_calibration_2805(spots_um, sphere_um)
    keep, diag = select_tracks_2805(tracks_um)
    print(f"  tracks kept    : {diag['n_kept']:,} ({100 * diag['n_kept']/diag['n_total']:.1f}%)")

    print("\nComputing per-step velocity...")
    spots_v_all = compute_per_step_velocity(spots)
    # Add vz (depth axis velocity) — needed for the YZ cross-section but not
    # produced by the 2508 module which only deals with in-plane (XY).
    if "vz" not in spots_v_all.columns:
        spots_v_all["vz"] = spots_v_all["dz"] / spots_v_all["dt_f"]
    n_dropped_speed = spots_v_all.attrs.get("n_dropped_speed", 0)
    spots_v = spots_v_all[spots_v_all["TRACK_ID"].isin(keep)]
    print(f"  per-step rows  : {len(spots_v):,}  "
          f"(dropped {n_dropped_speed:,} for instant > "
          f"{MAX_INSTANT_SPEED_UM_PER_FRAME} µm/frame)")

    arrow_scale = 40  # µm-per-frame per arrow length (in data units)
    arrow_lw = 1.3

    print(f"\nGenerating cross-section for frames {INTERN_FRAMES} "
          f"(window ±10 frames)...")
    fig = build_figure(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=INTERN_FRAMES, bin_y=BIN_Y_UM, bin_z=BIN_Z_UM,
        arrow_scale=arrow_scale, arrow_lw=arrow_lw, half_window=10,
    )
    with PdfPages(OUT_PDF) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    fig2 = build_figure(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=INTERN_FRAMES, bin_y=BIN_Y_UM, bin_z=BIN_Z_UM,
        arrow_scale=arrow_scale, arrow_lw=arrow_lw, half_window=10,
    )
    fig2.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    plt.close(fig2)

    print("Per-frame bin summary:")
    rows = build_summary_rows(spots_v, keep, sphere, INTERN_FRAMES, half_window=10)
    summary_df = pd.DataFrame(rows)
    print(summary_df.to_string(index=False))
    summary_df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved PDF to {OUT_PDF}")
    print(f"Saved PNG to {OUT_PNG}")
    print(f"Saved summary to {OUT_CSV}")


if __name__ == "__main__":
    main()