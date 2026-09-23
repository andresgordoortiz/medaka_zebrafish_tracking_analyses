"""Nuclei density at the margin (medaka 2508), single-frame snapshot.

Sibling to ``mk2508_surface_flow_maps_1h_back.py`` (which uses a 1h
backward window for velocity estimation). For density, we use a
SINGLE-FRAME snapshot at each key frame: averaging density over a 1h
window would smear the count by 60x and give nonsensical per-mm^2 numbers.

For each key frame, we:
  1. Take a single-frame snapshot of every nucleus at that FRAME.
  2. Localise the *margin* band by Y in data coordinates (animal pole is at
     negative Y after ``orient_embryo()``): **POSITION_Y in
     [MARGIN_Y_MIN_UM, MARGIN_Y_MAX_UM] = [-200, -100] um** (provided by
     the user from inspecting the vector maps).
  3. Restrict to the dorsal-ventral midline sagittal band (|Z - Z_MID|
     <= Z_BAND_UM) and the embryo's outer surface shell
     (SPHERICAL_DEPTH <= SURFACE_UM AND RADIAL_DIST >= R - SURFACE_UM),
     so the density reflects nuclei *at the surface only* (no ingressing
     cells).
  4. Bin the margin into thin X-tiles (DENSITY_BIN_UM = 20 um), count
     nuclei per tile, and normalise to nuclei per mm^2 (slab =
     20 um x 100 um x 30 um).
  5. Render a 1-row strip of density profiles (one panel per key frame)
     and a Panel B trend (mean margin density vs frame).

The margin is a ~100 um wide ring (Y direction). In data Y coordinates,
the embryo's animal pole is at y = -R = -840 um and the vegetal pole at
y = +R = +840 um. So the band y in [-200, -100] sits 100..200 um from
the animal pole -- the *embryonic margin* (where epiboly is active).

Outputs (in this directory):
  mk2508_margin_density_1h_back.pdf     -- page 1 = density strips, page 2 = trend
  mk2508_margin_density_1h_back.png     -- page 1 standalone PNG
  mk2508_margin_density_1h_back_summary.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

import mk2508_surface_flow_maps as m
from mk2508_surface_flow_maps import (
    KEY_FRAMES, FRAME_LABELS,
    VOXEL_XY_UM_PX, FRAME_INTERVAL_SEC, FRAMES_PER_MIN,
    MIN_TRACK_LEN, MAX_INSTANT_SPEED_UM_PER_FRAME,
    SURFACE_UM, Z_MID, Z_BAND_UM,
    apply_voxel_calibration, select_tracks,
    load_data,
)


# Constants kept for API symmetry with the flow-map sibling; ignored for density.
WIN_BACK_FRAMES = 60

# Biological labels for each key frame (matches the flow-map figure).
FRAME_LABELS_DENSITY = {
    116: "3 h before first contraction",
    160: "Start of internalisation",
    236: "2 h before first contraction",
    356: "Mid-epiboly / pre-contraction",
    476: "First contraction (end)",
}

# --- Margin localisation ----------------------------------------------------
# Provided by the user from visual inspection of the vector maps:
#   the margin ring sits at Y in [-200, -100] um in data coordinates
#   (animal pole = -840 um, vegetal pole = +840 um, embryo R = 840 um).
MARGIN_Y_MIN_UM = -200.0   # inclusive (farther from AP)
MARGIN_Y_MAX_UM = -100.0   # inclusive (closer to AP)
MARGIN_X_HALF_UM = 400.0   # +/-400 um wide X window around the dorsal-ventral
                           # midline (so we capture the full margin arc, not
                           # just its dorsal side)

# --- Binning ----------------------------------------------------------------
DENSITY_BIN_UM = 20.0      # um width of an X bin in the margin
MARGIN_WIDTH_UM = MARGIN_Y_MAX_UM - MARGIN_Y_MIN_UM  # = 100 um (the Y band)
DEPTH_BIN_UM    = SURFACE_UM                          # = 30 um (the surface cap)

HERE = Path(__file__).resolve().parent
IO_DIR = Path("results/medaka_mk2508")
STEM = "mk2508_margin_density_1h_back"
OUT_PDF = IO_DIR / f"{STEM}.pdf"
OUT_PNG = IO_DIR / f"{STEM}.png"
OUT_CSV = IO_DIR / f"{STEM}_summary.csv"


# =============================================================================
# Helpers
# =============================================================================
def select_window_back(spots: pd.DataFrame, frame: int,
                       valid_tracks: pd.Series,
                       win_back: int = WIN_BACK_FRAMES) -> pd.DataFrame:
    """Single-frame snapshot at *frame* (positional density, no time average).

    Density is a positional count: one row per nucleus per FRAME. Averaging
    over a 1h backward window (60 frames) would smear the count by 60x and
    give nonsensical per-mm^2 numbers. So we keep only rows with
    FRAME == key_frame (a single-frame snapshot). The ``win_back`` parameter
    is accepted for API symmetry with the flow-map sibling but is ignored.

    NOTE: We use the voxel-calibrated ``spots`` dataframe directly (NOT
    the per-step velocity table). The per-step speed filter (>3 um/frame)
    drops up to 50% of rows at the first contraction frame (476), which
    would artificially deflate the density. For positional counting, no
    speed filter is needed.
    """
    return spots[
        (spots["FRAME"] == frame)
        & spots["TRACK_ID"].isin(valid_tracks)
    ]


def select_margin(window: pd.DataFrame,
                  y_min_um: float = MARGIN_Y_MIN_UM,
                  y_max_um: float = MARGIN_Y_MAX_UM,
                  x_half_um: float = MARGIN_X_HALF_UM,
                  z_mid: float = Z_MID,
                  z_band_um: float = Z_BAND_UM,
                  surface_um: float = SURFACE_UM,
                  sphere_r_um: float | None = None) -> pd.DataFrame:
    """Pick nuclei in the margin band + sagittal slice + outer surface shell.

    Two cuts ensure we count **surface cells only** (no ingressing cells):

      - ``SPHERICAL_DEPTH <= surface_um``: cell is in the outermost cap
        layer (depth measured from the surface inward). Ingressing cells
        have larger SPHERICAL_DEPTH.
      - ``RADIAL_DIST >= sphere_r - surface_um``: cell is in the outermost
        radial shell. Cells that have moved radially inward (i.e., below
        the surface cap) have smaller RADIAL_DIST.

    The radial shell is a redundancy with SPHERICAL_DEPTH; both are kept
    because SPHERICAL_DEPTH can be negative (cells above the cap, where
    RADIAL_DIST > R). Combined they cleanly exclude ingressing cells.

    All cuts are in physical um (positions are calibrated before this call).
    ``y_max_um > y_min_um`` is required (data Y is negative at the animal
    pole, so [-200, -100] selects the band ~100-200 um away from the
    animal pole, toward the equator).
    """
    if y_max_um < y_min_um:
        y_min_um, y_max_um = y_max_um, y_min_um
    mask = (
        window["POSITION_Y"].between(y_min_um, y_max_um)
        & (window["POSITION_X"].abs() <= x_half_um)
        & ((window["POSITION_Z"] - z_mid).abs() <= z_band_um)
        & (window["SPHERICAL_DEPTH"] <= surface_um)
    )
    if sphere_r_um is not None:
        mask = mask & (window["RADIAL_DIST"] >= sphere_r_um - surface_um)
    return window[mask].copy()


def density_profile(margin: pd.DataFrame, bin_um: float = DENSITY_BIN_UM,
                     margin_width_um: float = MARGIN_WIDTH_UM):
    """Compute per-bin snapshot density along the margin's X axis.

    The ``margin`` dataframe contains one row per nucleus at the snapshot
    frame, restricted to the outermost surface shell
    (SPHERICAL_DEPTH <= 30 AND RADIAL_DIST >= R - 30).

    For each X bin:

        count           = #nuclei in the bin at this frame
        footprint_um2   = bin_um * margin_width_um         (X width * Y band)
        density_per_mm2 = count / footprint_um2 * 1e6      (nuclei / mm^2)

    The depth dimension is already accounted for by restricting rows to
    the outermost surface shell -- we are measuring surface *area*
    density, not volumetric density. The denominator is therefore the
    2D footprint area (X * Y), NOT the slab volume (X * Y * Z).

    Returns columns: x_bin, count, n_frames, mean_count_per_frame (=count),
    density_per_mm2, frac.
    """
    if len(margin) == 0:
        return pd.DataFrame(columns=["x_bin", "count", "n_frames",
                                     "mean_count_per_frame",
                                     "density_per_mm2", "frac"])

    df = margin.copy()
    df["x_bin"] = np.floor(df["POSITION_X"] / bin_um) * bin_um + bin_um / 2
    g = (df.groupby("x_bin", as_index=False)
           .size()
           .rename(columns={"size": "count"}))
    g["n_frames"] = df["FRAME"].nunique()
    g["mean_count_per_frame"] = g["count"].astype(float)  # single-frame snapshot
    total = float(g["count"].sum())
    g["frac"] = g["count"] / max(total, 1.0)
    footprint_um2 = bin_um * margin_width_um
    g["density_per_mm2"] = g["count"].astype(float) / footprint_um2 * 1e6
    return g.sort_values("x_bin").reset_index(drop=True)


# =============================================================================
# Plot
# =============================================================================
def draw_density_panel(ax, prof: pd.DataFrame, frame: int,
                       bin_um: float,
                       xlim_um: tuple[float, float]) -> None:
    """Bar plot of nuclei density vs X for one key frame."""
    if len(prof) == 0:
        ax.text(0.5, 0.5, "no nuclei in margin", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(f"Frame {frame}\n{FRAME_LABELS_DENSITY.get(frame, '')}",
                     fontsize=11, fontweight="bold")
        ax.set_xlim(*xlim_um)
        return

    ax.bar(prof["x_bin"], prof["density_per_mm2"],
           width=bin_um * 0.95, color="#4393C3", edgecolor="black",
           linewidth=0.4, alpha=0.85)

    ax.axvline(0, color="black", lw=0.6, alpha=0.5)
    ax.text(0.01, 0.97,
            f"n_obs = {int(prof['count'].sum()):,}\n"
            f"mean = {prof['density_per_mm2'].mean():.1f} /mm^2",
            transform=ax.transAxes, fontsize=8, va="top", ha="left",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))
    ax.set_xlim(*xlim_um)
    ax.set_xlabel("X (dorsal -> +X)", fontsize=10)
    ax.set_title(f"Frame {frame}\n{FRAME_LABELS_DENSITY.get(frame, '')}",
                 fontsize=11, fontweight="bold")


def _build_page1(*, density_by_frame: dict[int, pd.DataFrame],
                 xlim_um: tuple[float, float]) -> plt.Figure:
    """Page 1: density profile per key frame (one panel per frame)."""
    frames = list(density_by_frame.keys())
    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(5.0 * n, 5.2), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, f in zip(axes, frames):
        draw_density_panel(ax, density_by_frame[f], frame=f,
                           bin_um=DENSITY_BIN_UM,
                           xlim_um=xlim_um)

    axes[0].set_ylabel(f"Surface nuclei density at the margin\n"
                       f"(nuclei / mm^2, footprint = "
                       f"{DENSITY_BIN_UM:.0f} um X-bin x {MARGIN_WIDTH_UM:.0f} um Y-band)",
                       fontsize=10)

    fig.suptitle(
        "Medaka 2508 - Surface nuclei density at the embryonic margin "
        "(single-frame snapshot at each key frame)\n"
        f"Margin: Y in [{MARGIN_Y_MIN_UM:.0f}, {MARGIN_Y_MAX_UM:.0f}] um  |  "
        f"Surface shell: SPHERICAL_DEPTH <= {SURFACE_UM:.0f} um AND "
        f"RADIAL_DIST >= R - {SURFACE_UM:.0f} um  |  voxel {VOXEL_XY_UM_PX} um/px",
        fontsize=10.5, fontweight="bold", y=1.04,
    )
    fig.tight_layout()
    return fig


def _build_page2(stats: list[dict]) -> plt.Figure:
    """Page 2: trend across the 5 timepoints.

    Two stacked panels sharing the X axis:
      - Top: mean margin density vs frame index
      - Bottom: total nuclei observed in the margin band per frame
    """
    df = pd.DataFrame(stats).sort_values("frame").reset_index(drop=True)
    frames = df["frame"].to_numpy()
    density = df["mean_density_per_mm2"].to_numpy()
    counts = df["n_nuclei_obs"].to_numpy()

    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    fig.suptitle(
        "Medaka 2508 - Margin density trend across timepoints "
        "(single-frame snapshot at each key frame)\n"
        f"tracks >= {MIN_TRACK_LEN} spots | per-step instant <= "
        f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.0f} um/frame | "
        f"footprint = {DENSITY_BIN_UM:.0f} um X-bin x {MARGIN_WIDTH_UM:.0f} um Y-band "
        f"(surface shell: SPHERICAL_DEPTH <= {SURFACE_UM:.0f} um AND "
        f"RADIAL_DIST >= R - {SURFACE_UM:.0f} um)",
        fontsize=10.5, fontweight="bold", y=1.01,
    )

    # --- Top: mean density ---
    ax = axes[0]
    ax.plot(frames, density, "o-", color="#1F1F1F", lw=1.6, ms=10,
            markerfacecolor="#4393C3", markeredgecolor="black",
            markeredgewidth=0.6, zorder=3)
    for x, y in zip(frames, density):
        ax.annotate(f"{y:.0f}", xy=(x, y), xytext=(0, 8),
                    textcoords="offset points",
                    fontsize=9, ha="center", color="#2166AC", fontweight="bold")
    ax.set_ylabel(f"Mean margin density (nuclei / mm^2)\n"
                  f"footprint = {DENSITY_BIN_UM:.0f} um X-bin x {MARGIN_WIDTH_UM:.0f} um Y-band",
                  fontsize=10)
    ax.grid(alpha=0.25)

    # --- Bottom: total nuclei observed ---
    ax = axes[1]
    ax.bar(frames, counts, width=20, color="#4393C3", edgecolor="black",
           linewidth=0.6, alpha=0.85)
    for x, y in zip(frames, counts):
        ax.annotate(f"{int(y):,}", xy=(x, y), xytext=(0, 5),
                    textcoords="offset points",
                    fontsize=9, ha="center", color="#1F1F1F", fontweight="bold")
    ax.set_ylabel("Total nuclei observed\n(snapshot at each key frame)", fontsize=10)
    ax.set_xlabel("Frame index", fontsize=11)
    ax.set_xticks(frames)
    ax.set_xticklabels([str(f) for f in frames], fontsize=9)
    ax.grid(alpha=0.25, axis="y")

    # Per-frame biological labels along the bottom
    y_lo, y_hi = ax.get_ylim()
    for x, f in zip(frames, [int(f) for f in frames]):
        label = FRAME_LABELS_DENSITY.get(f, "")
        ax.text(x, y_lo - (y_hi - y_lo) * 0.08,
                label, fontsize=8, ha="center", va="top", color="grey", wrap=True)

    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return fig


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    io_dir = IO_DIR

    print(f"\n{'='*70}\n  MEDAKA 2508 MARGIN DENSITY (single-frame snapshot)\n{'='*70}\n")

    print("Loading data...")
    spots_px, tracks, sphere_px = load_data(io_dir)
    print(f"  oriented_spots : {len(spots_px):,} rows (raw, in pixels)")
    print(f"  filtered_tracks: {len(tracks):,} rows")
    print(f"  sphere center  : ({sphere_px['center_x']:.1f}, {sphere_px['center_y']:.1f}, "
          f"{sphere_px['center_z']:.1f}) px  R = {sphere_px['radius']:.1f} px")

    print(f"\nVoxel calibration: {VOXEL_XY_UM_PX} um/px (XY), "
          f"{FRAME_INTERVAL_SEC}s/frame -> {FRAMES_PER_MIN} frames/min")
    spots, sphere = apply_voxel_calibration(spots_px, sphere_px)
    print(f"  sphere R in um: {sphere['radius']:.1f} um")

    keep, diag = select_tracks(tracks)
    print("\nTrack filtering (comet / meteor removal):")
    print(f"  total tracks        : {diag['n_total']:,}")
    print(f"  with n_spots >= {MIN_TRACK_LEN:<3} : {diag['n_long']:,}")
    print(f"  tracks removed as comet: {diag['n_comet_removed']:,}")
    print(f"  tracks kept          : {diag['n_kept']:,}  "
          f"({100 * diag['n_kept'] / diag['n_total']:.1f}%)")

    print("\nFor density: using voxel-calibrated spots directly (no speed filter).")
    print("  The per-step speed filter (>3 um/frame) drops up to 50% of rows")
    print("  at the first contraction frame, which would artificially deflate")
    print("  density. For positional counting, we keep all rows from kept tracks.")
    spots_for_density = spots[spots["TRACK_ID"].isin(keep)].copy()

    # --- Margin localisation confirmation ---
    print(f"\nMargin localisation: Y in [{MARGIN_Y_MIN_UM:.0f}, "
          f"{MARGIN_Y_MAX_UM:.0f}] um (animal pole at -R = "
          f"{-sphere['radius']:.0f} um, vegetal pole at +R = +{sphere['radius']:.0f} um)")
    print(f"  X window: +/-{MARGIN_X_HALF_UM:.0f} um around the DV midline")
    print(f"  Z band:   +/-{Z_BAND_UM:.0f} um (sagittal slice)")
    print(f"  Surface shell (excludes ingressing cells):")
    print(f"    SPHERICAL_DEPTH <= {SURFACE_UM:.0f} um  AND  "
          f"RADIAL_DIST >= {sphere['radius']:.0f} - {SURFACE_UM:.0f} = "
          f"{sphere['radius'] - SURFACE_UM:.0f} um")

    # --- Per-frame density profiles ---
    print(f"\nComputing margin density for frames {KEY_FRAMES} "
          f"(single-frame snapshots, surface only)...")
    density_by_frame: dict[int, pd.DataFrame] = {}
    stats = []
    for f in KEY_FRAMES:
        win = select_window_back(spots_for_density, f, keep, WIN_BACK_FRAMES)
        margin = select_margin(win, sphere_r_um=sphere["radius"])
        prof = density_profile(margin)
        density_by_frame[f] = prof
        if len(prof) > 0:
            stat = {
                "frame": f,
                "label": FRAME_LABELS_DENSITY[f],
                "n_nuclei_obs": int(prof["count"].sum()),
                "n_frames": int(prof["n_frames"].max()),
                "n_bins": int(len(prof)),
                "n_nuclei_per_frame": float(prof["mean_count_per_frame"].sum()),
                "mean_density_per_mm2": float(prof["density_per_mm2"].mean()),
                "median_density_per_mm2": float(prof["density_per_mm2"].median()),
                "max_density_per_mm2": float(prof["density_per_mm2"].max()),
                "total_nuclei_per_mm2_band": float(prof["density_per_mm2"].sum()),
            }
        else:
            stat = {
                "frame": f, "label": FRAME_LABELS_DENSITY[f],
                "n_nuclei_obs": 0, "n_frames": 0, "n_bins": 0,
                "n_nuclei_per_frame": np.nan,
                "mean_density_per_mm2": np.nan, "median_density_per_mm2": np.nan,
                "max_density_per_mm2": np.nan, "total_nuclei_per_mm2_band": np.nan,
            }
        stats.append(stat)
        print(f"  frame {f:>3}: snapshot at frame {f}   "
              f"n_nuclei = {stat['n_nuclei_obs']:>5,}   "
              f"mean density = {stat['mean_density_per_mm2']:>7.0f} /mm^2   "
              f"({stat['n_bins']} X-bins)")

    summary_df = pd.DataFrame(stats)
    summary_df.to_csv(OUT_CSV, index=False)

    # --- Page 1: density strip ---
    R_um = sphere["radius"]
    xlim_um = (-R_um * 0.6, R_um * 0.6)

    page1 = _build_page1(density_by_frame=density_by_frame, xlim_um=xlim_um)

    # --- Page 2: trend ---
    page2 = _build_page2(stats)

    # --- Save ---
    with PdfPages(OUT_PDF) as pdf:
        pdf.savefig(page1, bbox_inches="tight")
        pdf.savefig(page2, bbox_inches="tight")
    plt.close(page1)
    plt.close(page2)

    # Standalone PNG of page 1
    page1b = _build_page1(density_by_frame=density_by_frame, xlim_um=xlim_um)
    page1b.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    plt.close(page1b)

    print(f"\nSaved 2-page PDF to {OUT_PDF}")
    print(f"Saved page-1 PNG to {OUT_PNG}")
    print(f"Saved summary to {OUT_CSV}")

    print("\nPer-frame summary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
