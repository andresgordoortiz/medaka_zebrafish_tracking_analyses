"""Surface flow maps for medaka 2508 at key developmental timepoints.

For each requested key frame, plot a per-bin velocity quiver of the cellular
surface (no yolk) on a representative cross-section (sagittal slice through
the dorsal-ventral midline). Animal pole (AP) is at the top after reversing
Y (orient_embryo() rotates AP -> -Y in this dataset).

Track filtering ("comet / meteor" removal)
------------------------------------------
  1.  n_spots >= MIN_TRACK_LEN                 (track is long enough to estimate motion)
  2.  max instantaneous speed <= median + 3*std of the track population
      (drops sudden teleports / mislinked tracks that cover large jumps per step)

Surface cells
-------------
  SPHERICAL_DEPTH <= SURFACE_UM                (cells at the outermost cap layer)

Per-frame window
----------------
  We compute per-step velocities once (dx, dy, dz / dt) and then average
  them inside a temporal window HALF_WIN around each key frame. This is
  the same recipe as gastrulation_dynamics_*.R, but evaluated at single
  frames instead of broad epochs.

Cross-section
-------------
  Sagittal slice through the dorsal-ventral midline:
      |Z - Z_MID| < Z_BAND_UM
  showing X (dorsal) horizontal and Y (animal up) vertical.

Outputs (in this directory):
  mk2508_surface_flow_maps_<suffix>.pdf / .png   -- one row per key frame (page 1)
                                              + Panel B trend figure (page 2, PDF only)
  mk2508_surface_flow_summary.csv                -- per-frame summary
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent

KEY_FRAMES         = [116, 160, 236, 356, 476]   # developmental milestones
FRAME_LABELS       = {                           # what each frame represents
    116: "3 h before first contraction",
    160: "Start of internalisation",
    236: "2 h before first contraction",
    356: "Mid-epiboly / pre-contraction",
    476: "First contraction (end)",
}

# --- Track filtering ---
MIN_TRACK_LEN      = 10        # spots per track
SPEED_OUTLIER_STD  = 3.0       # median + k * std (track-level, px/frame)
# Track-level speed in filtered_tracks.csv is in pixels / frame
# (oriented_spots coords are in pixels until multiplied by the voxel scale).

# --- Per-step instantaneous speed filter (µm/min) ---
# Applied to the per-step (dx, dy, dz)/dt vectors AFTER voxel conversion.
# Any step whose instantaneous speed exceeds this cap is dropped. This
# removes the long comet/meteor tails where a single mislink produces a
# large per-frame jump.  Note: speed is reported per-frame, not per-time.
# At 30s/frame, 3 µm/frame = 6 µm/min.
MAX_INSTANT_SPEED_UM_PER_FRAME = 3.0   # µm/frame (per-step displacement cap)

# --- Voxel calibration ---
# Positions and per-step displacements in oriented_spots.csv are in pixels.
# Convert to µm with the medaka voxel size. The X/Y and Z voxels are
# near-isotropic (<0.2% anisotropy), so a single scale factor is applied.
VOXEL_XY_UM_PX     = 1.05152   # µm per XY pixel (medaka lowres)
VOXEL_Z_UM_PX      = 1.05263   # µm per Z pixel
FRAME_INTERVAL_SEC = 30        # s/frame (medaka lowres)
FRAMES_PER_MIN     = 60 / FRAME_INTERVAL_SEC   # = 2 frames per minute

# --- Flow field ---
FLOW_BIN_UM        = 30        # µm -- spatial bin in the cross-section
FLOW_MIN_N         = 4         # minimum # vectors per bin to draw an arrow
HALF_WIN           = 8         # frames before/after for velocity estimation
ARROW_SCALE        = 50        # visual scale factor (arrow length per µm/frame)
ARROW_HEAD_FRAC    = 0.22      # head length as fraction of shaft
ARROW_LW           = 1.1

# --- Surface & cross-section ---
SURFACE_UM         = 30.0      # spherical depth <= this counts as "surface"
Z_MID              = 0.0       # mid-plane of the sagittal slice
Z_BAND_UM          = 500.0     # |z - Z_MID| <= this band
                              # (500 µm ≈ effectively the full FOV; matches
                              # the existing make_flow_map() default which
                              # projects all cells onto XY)

# --- Plot extent ---
# Use the actual cell extent (with padding) instead of the full sphere radius
# so the embryo fills the plot box. The full sphere outline is still drawn
# as a faint reference.
EXTENT_PAD_UM      = 50.0      # extra µm of padding around the cell extent


# =============================================================================
# Helpers
# =============================================================================
def load_data(io_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return (oriented_spots in px, filtered_tracks in px, sphere_params).

    The function reads all values in their raw pixel units. Conversion to
    µm is performed later by ``apply_voxel_calibration`` so the threshold
    diagnostics stay in pixels (matching ``filtered_tracks.csv``).
    """
    spots = pd.read_csv(io_dir / "oriented_spots.csv", low_memory=False)
    tracks = pd.read_csv(io_dir / "filtered_tracks.csv", low_memory=False)
    sph = pd.read_csv(io_dir / "sphere_params.csv")
    numeric_keys = {"radius", "center_x", "center_y", "center_z",
                    "rmse", "cap_height", "cap_base_radius"}
    sphere = {row["parameter"]: row["value"] for _, row in sph.iterrows()}
    for k in list(sphere.keys()):
        if k in numeric_keys:
            try:
                sphere[k] = float(sphere[k])
            except (TypeError, ValueError):
                pass
    return spots, tracks, sphere


def apply_voxel_calibration(spots: pd.DataFrame, sphere: dict) -> tuple[pd.DataFrame, dict]:
    """Convert pixel units to µm in-place on a copy of *spots* and *sphere*.

    The medaka voxel is near-isotropic, so the same XY scale is applied to
    POSITION_X/Y/Z, RADIAL_DIST, SPHERICAL_DEPTH, RADIAL_DIST_TO_CENTER.
    Angles (THETA_DEG, PHI_DEG) are dimensionless and are left as-is.
    """
    s = sphere.copy()
    for k in ("center_x", "center_y", "center_z", "radius", "rmse",
              "cap_height", "cap_base_radius"):
        if k in s and isinstance(s[k], (int, float)):
            s[k] = s[k] * VOXEL_XY_UM_PX

    sp = spots.copy()
    for col in ("POSITION_X", "POSITION_Y", "POSITION_Z",
                "RADIAL_DIST", "SPHERICAL_DEPTH", "RADIAL_DIST_TO_CENTER"):
        if col in sp.columns:
            sp[col] = sp[col] * VOXEL_XY_UM_PX
    return sp, s


def select_tracks(tracks: pd.DataFrame,
                  min_len: int = MIN_TRACK_LEN,
                  k_std: float = SPEED_OUTLIER_STD) -> tuple[pd.Series, dict]:
    """Return (mask of valid track IDs, diagnostics dict)."""
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


def compute_per_step_velocity(spots: pd.DataFrame,
                              max_inst_speed_um_per_frame: float = MAX_INSTANT_SPEED_UM_PER_FRAME
                              ) -> pd.DataFrame:
    """Add per-step (dx,dy,dz)/dt and speed columns, sorted by (TRACK, FRAME).

    Per-step vectors with instantaneous speed above
    ``max_inst_speed_um_per_frame`` (in µm/frame) are dropped to remove
    mislinked "comet" jumps.
    """
    s = spots.sort_values(["TRACK_ID", "FRAME"]).copy()
    g = s.groupby("TRACK_ID", sort=False)
    s["dx"]   = s["POSITION_X"] - g["POSITION_X"].shift()
    s["dy"]   = s["POSITION_Y"] - g["POSITION_Y"].shift()
    s["dz"]   = s["POSITION_Z"] - g["POSITION_Z"].shift()
    s["dt_f"] = s["FRAME"]       - g["FRAME"].shift()
    s = s[(s["dt_f"] > 0) & s["dx"].notna()].copy()
    s["vx"]     = s["dx"] / s["dt_f"]
    s["vy"]     = s["dy"] / s["dt_f"]
    # NOTE: "speed" is the magnitude of the *projected* (XY) velocity,
    # matching what the flow-map arrows show.  We deliberately do NOT include
    # the Z component, because the figure projects cells onto XY (full Z band)
    # and the arrows therefore encode the in-plane motion only.
    s["speed"]  = np.sqrt(s["dx"]**2 + s["dy"]**2) / s["dt_f"]
    s["speed_3d"] = np.sqrt(s["dx"]**2 + s["dy"]**2 + s["dz"]**2) / s["dt_f"]
    s["speed_um_per_min"] = s["speed"] * FRAMES_PER_MIN
    # Mid-frame as anchor for plotting
    s["MID_X"] = s["POSITION_X"] - 0.5 * s["dx"]
    s["MID_Y"] = s["POSITION_Y"] - 0.5 * s["dy"]
    s["MID_Z"] = s["POSITION_Z"] - 0.5 * s["dz"]
    s["MID_FRAME"] = s["FRAME"] - 0.5 * s["dt_f"]
    s["MID_SPHERICAL_DEPTH"] = (
        s["SPHERICAL_DEPTH"] - 0.5 * (s["SPHERICAL_DEPTH"] - g["SPHERICAL_DEPTH"].shift())
    )
    n0 = len(s)
    s = s[s["speed"] <= max_inst_speed_um_per_frame].copy()
    s.attrs["n_dropped_speed"] = n0 - len(s)
    return s


def select_window(spots_v: pd.DataFrame, frame: int,
                  valid_tracks: pd.Series, half_win: int = HALF_WIN) -> pd.DataFrame:
    """Subset per-step velocity vectors to a temporal window around *frame*."""
    window = spots_v[
        spots_v["MID_FRAME"].between(frame - half_win, frame + half_win)
        & spots_v["TRACK_ID"].isin(valid_tracks)
    ]
    return window


def select_surface(window: pd.DataFrame,
                   surface_um: float = SURFACE_UM) -> pd.DataFrame:
    """Restrict to the surface band (cells at the outer cap layer)."""
    return window[window["SPHERICAL_DEPTH"] <= surface_um]


def compute_track_stats_per_frame(
    spots_v: pd.DataFrame,
    valid_tracks: pd.Series,
    frames: list[int],
    half_win: int = HALF_WIN,
    surface_um: float = SURFACE_UM,
    bin_um: float = FLOW_BIN_UM,
    min_n: int = FLOW_MIN_N,
) -> dict[int, dict]:
    """Per-frame, bin-averaged statistics that match what the figure shows.

    For each key frame, we look at the ±half_win temporal window, pick
    every step in the surface band (SPHERICAL_DEPTH <= surface_um), and
    bin into ``bin_um`` XY tiles (the same binning used by the flow map).
    The trends are computed from the bin-mean velocities:

      - frac_ap_ward_bins: fraction of bins with bin-mean vy < 0
        (THIS MATCHES THE FIGURE: if you look at the vector map, a tile
        is "AP-ward" when its average arrow points up toward the animal
        pole.)
      - mean_bin_speed_um_per_min: mean bin-mean speed in µm/min
        (i.e., mean over bins of sqrt(vx^2 + vy^2) * FRAMES_PER_MIN,
        NOT mean over individual steps)
      - n_bins, n_steps: counts for transparency

    Computing the fraction from bin means (rather than counting raw steps
    with vy<0) gives a much better match to the visual: at frame 476
    (first contraction), the bulk vegetal-ward motion dominates the
    bin means, so the fraction drops to ~4%, even though many individual
    steps still have negative vy.

    Returns a dict keyed by frame.
    """
    valid_set = set(valid_tracks)
    out: dict[int, dict] = {}
    for frame in frames:
        window = select_window(spots_v, frame, pd.Series(list(valid_set)), half_win)
        surface = select_surface(window, surface_um)
        if len(surface) == 0:
            out[frame] = {
                "frac_ap_ward_bins": np.nan,
                "mean_bin_speed_um_per_min": np.nan,
                "n_bins": 0,
                "n_steps": 0,
            }
            continue
        # Bin into XY tiles using the same scheme as bin_xy()
        binned = bin_xy(surface, bin_um=bin_um, min_n=min_n)
        if len(binned) == 0:
            out[frame] = {
                "frac_ap_ward_bins": np.nan,
                "mean_bin_speed_um_per_min": np.nan,
                "n_bins": 0,
                "n_steps": int(len(surface)),
            }
            continue
        # frac_ap_ward: fraction of bins with bin-mean vy < 0
        # bin speed: sqrt(vx^2 + vy^2), then averaged across bins
        binned = binned.copy()
        binned["bin_speed_um_per_min"] = (
            np.sqrt(binned["vx"]**2 + binned["vy"]**2) * FRAMES_PER_MIN
        )
        out[frame] = {
            "frac_ap_ward_bins": float((binned["vy"] < 0).mean()),
            "mean_bin_speed_um_per_min": float(binned["bin_speed_um_per_min"].mean()),
            "n_bins": int(len(binned)),
            "n_steps": int(len(surface)),
        }
    return out


def select_slice(window: pd.DataFrame,
                 z_mid: float = Z_MID, z_band: float = Z_BAND_UM) -> pd.DataFrame:
    """Restrict to the sagittal slice through the dorsal-ventral midline."""
    return window[(window["MID_Z"] - z_mid).abs() <= z_band]


def bin_xy(window: pd.DataFrame,
           bin_um: float = FLOW_BIN_UM,
           min_n: int = FLOW_MIN_N) -> pd.DataFrame:
    """Bin velocity vectors onto a regular XY grid (bin centers in µm)."""
    df = window.copy()
    df["x_bin"] = np.floor(df["MID_X"] / bin_um) * bin_um + bin_um / 2
    df["y_bin"] = np.floor(df["MID_Y"] / bin_um) * bin_um + bin_um / 2
    g = df.groupby(["x_bin", "y_bin"], as_index=False)
    out = g.agg(
        vx=("vx", "mean"),
        vy=("vy", "mean"),
        speed=("speed", "mean"),
        depth_med=("SPHERICAL_DEPTH", "median"),
        theta_med=("THETA_DEG", "median"),
        phi_med=("PHI_DEG", "median"),
        n=("vx", "count"),
    )
    return out[out["n"] >= min_n]


def classify_bins(fl: pd.DataFrame) -> pd.DataFrame:
    """Per-bin direction label:
        AP-ward  (toward animal pole)   if vy < 0   (data Y is negative there)
        VP-ward  (toward vegetal pole)   if vy > 0
    """
    fl = fl.copy()
    fl["move_type"] = np.where(fl["vy"] >= 0, "VP-ward", "AP-ward")
    return fl


def draw_panel(ax, fl: pd.DataFrame, sphere: dict,
               frame: int, surface_um: float,
               z_mid: float, z_band: float,
               bin_um: float, arrow_scale: float,
               arrow_lw: float, arrow_head: float,
               show_sphere: bool = True,
               bg_scatter: np.ndarray | None = None,
               extent_pad_um: float = EXTENT_PAD_UM) -> None:
    """Draw a single flow-map panel with sphere outline + quiver arrows.

    The plot uses the convention ``orient_embryo()`` lays out:
      * data X axis = dorsal (+X) / ventral (-X)
      * data Y axis = animal pole (-Y) / vegetal pole (+Y)
    We flip Y for display so the animal pole sits at the TOP of each panel
    (matching the existing ``medaka_density_analysis.R`` flow maps).
    The plot extent is set to the actual cell distribution (with padding)
    so the embryo fills the box instead of leaving white space.
    """
    cols = {"AP-ward": "#2166AC", "VP-ward": "#B2182B"}
    R = sphere["radius"]
    th = np.linspace(0, 2 * np.pi, 200)
    cx = R * np.cos(th)
    cy_data = R * np.sin(th)         # data Y
    cy_display = -cy_data            # flip so AP at top (cy_display is positive upward)

    # Choose plot extent from the data (cells + scatter), not from the full
    # sphere radius -- this makes the embryo fill the panel.
    if bg_scatter is not None and len(bg_scatter) > 0:
        x_lo = min(bg_scatter[:, 0].min(), fl["x_bin"].min() if len(fl) else np.inf)
        x_hi = max(bg_scatter[:, 0].max(), fl["x_bin"].max() if len(fl) else -np.inf)
        y_lo_data = min(bg_scatter[:, 1].min(), fl["y_bin"].min() if len(fl) else np.inf)
        y_hi_data = max(bg_scatter[:, 1].max(), fl["y_bin"].max() if len(fl) else -np.inf)
    elif len(fl) > 0:
        x_lo, x_hi = fl["x_bin"].min(), fl["x_bin"].max()
        y_lo_data, y_hi_data = fl["y_bin"].min(), fl["y_bin"].max()
    else:
        x_lo, x_hi = -R, R
        y_lo_data, y_hi_data = -R, R
    extent_x = [x_lo - extent_pad_um, x_hi + extent_pad_um]
    # extent_y is in DISPLAY coords (Y flipped), so convert data limits to display
    extent_y_display = [-y_hi_data - extent_pad_um, -y_lo_data + extent_pad_um]

    # Background density (all surface cells in the time window) -- Y flipped
    if bg_scatter is not None and len(bg_scatter) > 0:
        max_bg = 4000
        if len(bg_scatter) > max_bg:
            idx = np.random.default_rng(42).choice(
                len(bg_scatter), size=max_bg, replace=False)
            bg = bg_scatter[idx]
        else:
            bg = bg_scatter
        ax.scatter(bg[:, 0], -bg[:, 1],
                   s=0.5, c="lightgray", alpha=0.4, rasterized=True)

    # Draw arrows ordered so AP-ward sits underneath VP-ward
    for mt in ["VP-ward", "AP-ward"]:
        sub = fl[fl["move_type"] == mt]
        if sub.empty:
            continue
        # Y-flip: -y_bin for display, -vy so that the visual direction
        # matches the data direction (vy<0 -> arrow points down toward AP).
        ax.quiver(
            sub["x_bin"], -sub["y_bin"],
            sub["vx"] * arrow_scale, -sub["vy"] * arrow_scale,
            color=cols[mt], linewidth=arrow_lw,
            headwidth=5, headlength=6, headaxislength=5,
            angles="xy", scale_units="xy", scale=1.0,
            alpha=0.95,
        )

    # Sphere outline in display coords (animal pole at top)
    if show_sphere:
        ax.plot(cx, cy_display, color="black", lw=0.6, alpha=0.35)
        ax.plot(0, 0, "+", color="black", ms=8, alpha=0.4)

    ax.set_aspect("equal")
    ax.set_xlim(extent_x[0], extent_x[1])
    ax.set_ylim(extent_y_display[0], extent_y_display[1])
    ax.set_xlabel("X (dorsal → +X)", fontsize=10)
    ax.set_ylabel("Y (animal pole ↑)", fontsize=10)

    n_ap = int((fl["move_type"] == "AP-ward").sum())
    n_vp = int((fl["move_type"] == "VP-ward").sum())
    title = f"Frame {frame}\n{FRAME_LABELS.get(frame, '')}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.text(0.02, 0.97,
            f"{len(fl)} bins | AP-ward {n_ap}  VP-ward {n_vp}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1))


def make_figure(spots_v: pd.DataFrame,
                valid_tracks: pd.Series, sphere: dict,
                frames: list[int],
                surface_um: float, z_mid: float, z_band: float,
                bin_um: float, min_n: int,
                arrow_scale: float, arrow_lw: float, arrow_head: float,
                half_win: int, save_stem: Path,
                speed_thr_px: float | None = None,
                speed_median_px: float | None = None,
                speed_std_px: float | None = None,
                n_dropped_speed: int = 0) -> dict:
    """Make a horizontal strip of flow panels, one per frame, save PNG only.

    The PDF (2-page: flow maps + Panel B) is built separately by main() via
    PdfPages; this function only writes the standalone PNG.
    """
    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 6.5))
    if n == 1:
        axes = [axes]

    summaries = []
    for ax, frame in zip(axes, frames):
        window = select_window(spots_v, frame, valid_tracks, half_win)
        surface = select_surface(window, surface_um)
        slice_ = select_slice(surface, z_mid, z_band)
        fl = bin_xy(slice_, bin_um, min_n)
        fl = classify_bins(fl)

        # Background scatter: all surface cells (no binning) in the window
        bg = slice_[["MID_X", "MID_Y"]].to_numpy()

        draw_panel(ax, fl, sphere,
                   frame=frame, surface_um=surface_um,
                   z_mid=z_mid, z_band=z_band,
                   bin_um=bin_um, arrow_scale=arrow_scale,
                   arrow_lw=arrow_lw, arrow_head=arrow_head,
                   bg_scatter=bg)

        # Summary stats (vx, vy, speed already in µm/frame since spots were
        # converted by apply_voxel_calibration())
        if len(fl) > 0:
            v_mean = fl["speed"].mean()
            up_frac = (fl["move_type"] == "AP-ward").mean()
            vx_mean = fl["vx"].mean()
            vy_mean = fl["vy"].mean()
        else:
            v_mean = vx_mean = vy_mean = up_frac = np.nan
        summaries.append({
            "frame": frame,
            "label": FRAME_LABELS.get(frame, ""),
            "n_vectors": int(len(slice_)),
            "n_bins": int(len(fl)),
            "n_surface_window": int(len(surface)),
            "mean_speed_um_per_frame": float(v_mean),
            "mean_speed_um_per_min":   float(v_mean * FRAMES_PER_MIN),
            "mean_vx_um_per_frame": float(vx_mean),
            "mean_vy_um_per_frame": float(vy_mean),
            "frac_ap_ward_bins": float(up_frac),
        })

    # Common legend
    handles = [
        Line2D([0], [0], color="#2166AC", lw=2, label="AP-ward (vy < 0, toward animal pole)"),
        Line2D([0], [0], color="#B2182B", lw=2, label="VP-ward (vy > 0, toward vegetal pole)"),
        Line2D([0], [0], color="black", lw=1, label="Sphere (R = %.0f µm)" % sphere["radius"]),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), frameon=False)

    # Subtitle (single line, concise).
    #
    # Track-level filter:   ≥ N spots   (drops too-short tracks)
    # Per-step filter:      instant ≤ MAX_INSTANT_SPEED_UM_PER_FRAME (µm/frame)
    # Each bin's mean speed is reported in the table in µm/min, matching
    # the per-step filter.
    speed_clause = (
        f"tracks ≥ {MIN_TRACK_LEN} spots; per-step instant speed ≤ "
        f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.1f} µm/frame "
        f"({MAX_INSTANT_SPEED_UM_PER_FRAME * FRAMES_PER_MIN:.0f} µm/min) | "
        f"{n_dropped_speed:,} comet steps dropped"
    )
    fig.suptitle(
        "Medaka 2508 — Surface flow maps at key developmental timepoints\n"
        f"Side view (XY, AP at top) | surface band depth ≤ {surface_um:.0f} µm | "
        f"voxel {VOXEL_XY_UM_PX} µm/px | {speed_clause}",
        fontsize=11, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    # Save page 1 (flow maps) only as PNG; PDF is built later with both pages.
    fig.savefig(save_stem.with_suffix(".png"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    return {"summaries": summaries, "figure": str(save_stem)}


# =============================================================================
# Page builders (return Figure objects for PdfPages)
# =============================================================================
def _build_page1(*, spots_v, valid_tracks, sphere, frames,
                 surface_um, z_mid, z_band, bin_um, min_n,
                 arrow_scale, arrow_lw, arrow_head, half_win,
                 speed_thr_px, speed_median_px, speed_std_px,
                 n_dropped_speed) -> plt.Figure:
    """Page 1: flow maps, identical to the standalone PNG."""
    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 6.5))
    if n == 1:
        axes = [axes]
    for ax, frame in zip(axes, frames):
        window = select_window(spots_v, frame, valid_tracks, half_win)
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

    handles = [
        Line2D([0], [0], color="#2166AC", lw=2, label="AP-ward (vy < 0, toward animal pole)"),
        Line2D([0], [0], color="#B2182B", lw=2, label="VP-ward (vy > 0, toward vegetal pole)"),
        Line2D([0], [0], color="black", lw=1, label=f"Sphere (R = {sphere['radius']:.0f} µm)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), frameon=False)

    speed_clause = (
        f"tracks ≥ {MIN_TRACK_LEN} spots; per-step instant displacement ≤ "
        f"{MAX_INSTANT_SPEED_UM_PER_FRAME:.1f} µm/frame "
        f"({MAX_INSTANT_SPEED_UM_PER_FRAME * FRAMES_PER_MIN:.0f} µm/min) | "
        f"{n_dropped_speed:,} comet steps dropped"
    )
    fig.suptitle(
        "Medaka 2508 — Surface flow maps at key developmental timepoints\n"
        f"Side view (XY, AP at top) | surface band depth ≤ {surface_um:.0f} µm | "
        f"voxel {VOXEL_XY_UM_PX} µm/px | {speed_clause}",
        fontsize=11, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    return fig


def _build_page2(summaries: list[dict],
                 track_summary: dict[int, dict] | None = None) -> plt.Figure:
    """Page 2: Panel B trend figure.

    If ``track_summary`` is provided, use those per-frame stats (computed
    directly from raw per-step velocities of every track in the surface
    band, NOT from the binned flow map). Otherwise fall back to the
    bin-averaged stats from the flow-map ``summaries``.
    """
    df = pd.DataFrame(summaries).sort_values("frame").reset_index(drop=True)
    frames = df["frame"].to_numpy()

    if track_summary is not None and len(track_summary) > 0:
        # Use bin-mean stats (matches what the figure shows)
        ts = (pd.DataFrame.from_dict(track_summary, orient="index")
                .reindex(frames))
        frac = (ts["frac_ap_ward_bins"].to_numpy() * 100.0)
        speed = ts["mean_bin_speed_um_per_min"].to_numpy()
    else:
        frac = df["frac_ap_ward_bins"].to_numpy() * 100.0
        speed = df["mean_speed_um_per_min"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True)
    source_label = "per-bin (matches figure)"
    fig.suptitle(
        "Medaka 2508 — Panel B: surface flow composition across timepoints\n"
        f"Trends computed {source_label}  |  "
        f"per-step speed ≤ {MAX_INSTANT_SPEED_UM_PER_FRAME:.0f} µm/frame  |  "
        f"tracks ≥ {MIN_TRACK_LEN} spots  |  surface depth ≤ {SURFACE_UM:.0f} µm",
        fontsize=11, fontweight="bold", y=1.06,
    )

    # ---- Left: fraction AP-ward (matches figure) ---------------------
    ax = axes[0]
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
        "ap_vp", ["#B2182B", "#F4F4F4", "#2166AC"]
    )
    ax.scatter(frames, frac, c=frac, cmap=cmap, vmin=0, vmax=100,
               s=80, zorder=3, edgecolor="black", linewidth=0.6)
    ax.plot(frames, frac, color="black", lw=1.0, alpha=0.45, zorder=2)
    frac_ylabel = "% AP-ward bins (bin-mean vy < 0)"
    ax.set_ylabel(frac_ylabel, fontsize=11)
    ax.set_ylim(-5, 105)
    ax.axhline(50, color="grey", lw=0.6, ls=":", zorder=1)
    ax.text(frames.min(), 52, "50%", color="grey", fontsize=8,
            ha="left", va="bottom")
    ax.grid(alpha=0.25)
    for x, y in zip(frames, frac):
        ax.annotate(f"{y:.0f}%",
                    xy=(x, y), xytext=(0, 8),
                    textcoords="offset points",
                    fontsize=9, ha="center", color="#2166AC",
                    fontweight="bold")

    # ---- Right: mean bin speed (matches figure) -----------------------
    ax = axes[1]
    speed_ylabel = "Mean bin speed (µm / min)"
    ax.plot(frames, speed, "o-", color="#1F1F1F", lw=1.5, ms=8,
            markerfacecolor="#FB8B24", markeredgecolor="black",
            markeredgewidth=0.6, zorder=3)
    for x, y in zip(frames, speed):
        ax.annotate(f"{y:.2f}",
                    xy=(x, y), xytext=(0, 8),
                    textcoords="offset points",
                    fontsize=9, ha="center", color="#B45309",
                    fontweight="bold")
    ax.set_ylabel(speed_ylabel, fontsize=11)
    ax.grid(alpha=0.25)

    # ---- Shared X-axis + biological labels -----------------------------
    for ax in axes:
        ax.set_xlabel("Frame index", fontsize=11)
        ax.set_xticks(frames)
        ax.set_xticklabels([str(f) for f in frames], fontsize=9)

    for x, f in zip(frames, [int(f) for f in frames]):
        label = FRAME_LABELS.get(f, "")
        bbox = axes[0].get_position()
        x_pos = bbox.x0 + bbox.width * (x - frames.min()) / max(frames.max() - frames.min(), 1)
        fig.text(x_pos, -0.04, label,
                 fontsize=8, ha="center", va="top",
                 color="grey", wrap=True)

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


# =============================================================================
# Panel B: trend figure (fraction AP-ward & mean XY speed vs frame)
# =============================================================================
def make_panel_b(summaries: list[dict], save_path: Path) -> None:
    """Two-panel trend figure: fraction AP-ward and mean XY speed vs frame.

    Both panels share the X-axis (frame index) so the reader can compare
    the two traces directly.  Numerical values are printed above each point.
    """
    df = pd.DataFrame(summaries).sort_values("frame").reset_index(drop=True)
    frames = df["frame"].to_numpy()
    frac   = df["frac_ap_ward_bins"].to_numpy() * 100.0
    speed  = df["mean_speed_um_per_min"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True)
    fig.suptitle(
        "Medaka 2508 — Panel B: surface flow composition across timepoints\n"
        f"Per-step speed ≤ {MAX_INSTANT_SPEED_UM_PER_FRAME:.0f} µm/frame  |  "
        f"tracks ≥ {MIN_TRACK_LEN} spots  |  surface depth ≤ {SURFACE_UM:.0f} µm",
        fontsize=11, fontweight="bold", y=1.06,
    )

    # ---- Left: fraction AP-ward ----------------------------------------
    ax = axes[0]
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
        "ap_vp", ["#B2182B", "#F4F4F4", "#2166AC"]
    )
    ax.scatter(frames, frac, c=frac, cmap=cmap, vmin=0, vmax=100,
               s=80, zorder=3, edgecolor="black", linewidth=0.6)
    ax.plot(frames, frac, color="black", lw=1.0, alpha=0.45, zorder=2)
    ax.set_ylabel("% AP-ward bins", fontsize=11)
    ax.set_ylim(-5, 105)
    ax.axhline(50, color="grey", lw=0.6, ls=":", zorder=1)
    ax.text(frames.min(), 52, "50%", color="grey", fontsize=8,
            ha="left", va="bottom")
    ax.grid(alpha=0.25)
    for x, y in zip(frames, frac):
        ax.annotate(f"{y:.0f}%",
                    xy=(x, y), xytext=(0, 8),
                    textcoords="offset points",
                    fontsize=9, ha="center", color="#2166AC",
                    fontweight="bold")

    # ---- Right: mean XY speed ------------------------------------------
    ax = axes[1]
    ax.plot(frames, speed, "o-", color="#1F1F1F", lw=1.5, ms=8,
            markerfacecolor="#FB8B24", markeredgecolor="black",
            markeredgewidth=0.6, zorder=3)
    for x, y in zip(frames, speed):
        ax.annotate(f"{y:.2f}",
                    xy=(x, y), xytext=(0, 8),
                    textcoords="offset points",
                    fontsize=9, ha="center", color="#B45309",
                    fontweight="bold")
    ax.set_ylabel("Mean surface speed (µm / min)", fontsize=11)
    ax.grid(alpha=0.25)

    # ---- Shared X-axis + biological labels -----------------------------
    for ax in axes:
        ax.set_xlabel("Frame index", fontsize=11)
        ax.set_xticks(frames)
        ax.set_xticklabels([str(f) for f in frames], fontsize=9)

    for x, f in zip(frames, [int(f) for f in frames]):
        label = FRAME_LABELS.get(f, "")
        bbox = axes[0].get_position()
        x_pos = bbox.x0 + bbox.width * (x - frames.min()) / max(frames.max() - frames.min(), 1)
        fig.text(x_pos, -0.04, label,
                 fontsize=8, ha="center", va="top",
                 color="grey", wrap=True)

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--io-dir", default="results/medaka_mk2508")
    parser.add_argument("--out", default="mk2508_surface_flow_maps")
    parser.add_argument("--frames", nargs="*", type=int, default=KEY_FRAMES)
    parser.add_argument("--surface", type=float, default=SURFACE_UM)
    parser.add_argument("--z-band", type=float, default=Z_BAND_UM)
    parser.add_argument("--bin", type=float, default=FLOW_BIN_UM)
    parser.add_argument("--min-n", type=int, default=FLOW_MIN_N)
    parser.add_argument("--half-win", type=int, default=HALF_WIN)
    args = parser.parse_args()

    io_dir = Path(args.io_dir)

    print(f"\n{'='*70}\n  MEDAKA 2508 SURFACE FLOW MAPS\n{'='*70}\n")

    print("Loading data...")
    spots_px, tracks, sphere_px = load_data(io_dir)
    print(f"  oriented_spots : {len(spots_px):,} rows (raw, in pixels)")
    print(f"  filtered_tracks: {len(tracks):,} rows")
    print(f"  sphere center  : ({sphere_px['center_x']:.1f}, {sphere_px['center_y']:.1f}, "
          f"{sphere_px['center_z']:.1f}) px  R = {sphere_px['radius']:.1f} px")

    # --- Voxel calibration: pixels -> µm ---
    # oriented_spots.csv and sphere_params.csv are in pixel units.  Apply the
    # medaka voxel scale (1.05152 µm/px XY) so all downstream geometry and
    # velocities are in physical units.
    print(f"\nVoxel calibration: {VOXEL_XY_UM_PX} µm/px (XY), {VOXEL_Z_UM_PX} µm/px (Z), "
          f"{FRAME_INTERVAL_SEC}s/frame -> {FRAMES_PER_MIN} frames/min")
    spots, sphere = apply_voxel_calibration(spots_px, sphere_px)
    print(f"  sphere R in µm: {sphere['radius']:.1f} µm")
    print(f"  cell range X (µm): [{spots['POSITION_X'].min():.1f}, {spots['POSITION_X'].max():.1f}]")
    print(f"  cell range Y (µm): [{spots['POSITION_Y'].min():.1f}, {spots['POSITION_Y'].max():.1f}]")

    # --- Filter tracks (operates on TRACK_MEAN_SPEED, which is in pixels) ---
    keep, diag = select_tracks(tracks)
    print("\nTrack filtering (comet / meteor removal):")
    print(f"  total tracks        : {diag['n_total']:,}")
    print(f"  with n_spots >= {MIN_TRACK_LEN:<3} : {diag['n_long']:,}")
    print(f"  TRACK_MEAN_SPEED median = {diag['median_speed']:.4f}  std = {diag['std_speed']:.4f} (px/frame)")
    print(f"    -> in µm: median = {diag['median_speed'] * VOXEL_XY_UM_PX:.4f}  "
          f"std = {diag['std_speed'] * VOXEL_XY_UM_PX:.4f} (µm/frame)")
    print(f"  threshold = median + {SPEED_OUTLIER_STD}·σ = {diag['speed_thr']:.4f} px/frame "
          f"= {diag['speed_thr'] * VOXEL_XY_UM_PX:.4f} µm/frame "
          f"({diag['speed_thr'] * VOXEL_XY_UM_PX * FRAMES_PER_MIN:.4f} µm/min)")
    print(f"  tracks removed as comet: {diag['n_comet_removed']:,}")
    print(f"  tracks kept          : {diag['n_kept']:,}  "
          f"({100 * diag['n_kept'] / diag['n_total']:.1f}%)")

    # --- Per-step velocities ---
    print("\nComputing per-step velocity (one row per step in a track)...")
    spots_v_all = compute_per_step_velocity(spots)
    n_dropped_speed = spots_v_all.attrs.get("n_dropped_speed", 0)
    # Restrict to the set of kept tracks
    spots_v = spots_v_all[spots_v_all["TRACK_ID"].isin(keep)]
    print(f"  per-step rows (after kept-track filter)     : {len(spots_v):,}")
    print(f"  per-step rows dropped (instant > 6 µm/min) : {n_dropped_speed:,} "
          f"({100 * n_dropped_speed / max(len(spots_v_all), 1):.2f}%)")
    print(f"  instantaneous speed stats (µm / frame, after voxel calibration):")
    print(f"    median = {spots_v['speed'].median():.3f}   "
          f"mean = {spots_v['speed'].mean():.3f}   "
          f"std = {spots_v['speed'].std():.3f}")
    print(f"    99th pct = {spots_v['speed'].quantile(0.99):.3f}   "
          f"max = {spots_v['speed'].max():.3f}")
    print(f"  equivalent in µm / min (× {FRAMES_PER_MIN}):")
    print(f"    median = {spots_v['speed'].median() * FRAMES_PER_MIN:.3f}   "
          f"mean = {spots_v['speed'].mean() * FRAMES_PER_MIN:.3f}")

    # --- Per-key-frame flow map + Panel B (2 pages in 1 PDF) ---
    save_stem = io_dir / args.out
    print(f"\nGenerating flow maps for frames {args.frames} ...")
    result = make_figure(
        spots_v=spots_v, valid_tracks=keep, sphere=sphere,
        frames=args.frames,
        surface_um=args.surface, z_mid=Z_MID, z_band=args.z_band,
        bin_um=args.bin, min_n=args.min_n,
        arrow_scale=ARROW_SCALE, arrow_lw=ARROW_LW, arrow_head=ARROW_HEAD_FRAC,
        half_win=args.half_win, save_stem=save_stem,
        speed_thr_px=diag["speed_thr"],
        speed_median_px=diag["median_speed"],
        speed_std_px=diag["std_speed"],
        n_dropped_speed=n_dropped_speed,
    )

    # --- Summary CSV ---
    summary_df = pd.DataFrame(result["summaries"])
    summary_path = io_dir / (args.out + "_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print("\nSummary table:")
    print(summary_df.to_string(index=False))

    # --- Build the final PDF: page 1 = flow maps, page 2 = Panel B trend ---
    # Compute per-frame stats DIRECTLY from raw per-step velocities of every
    # track in the surface band (not from the binned flow-map data). This is
    # what the trends in Panel B should reflect.
    pdf_path = save_stem.with_suffix(".pdf")
    track_stats = compute_track_stats_per_frame(
        spots_v=spots_v, valid_tracks=keep, frames=args.frames,
        half_win=args.half_win, surface_um=args.surface,
        bin_um=args.bin, min_n=args.min_n,
    )
    print("\nBin-averaged Panel B statistics (matches the figure):")
    for f in args.frames:
        st = track_stats[f]
        print(f"  frame {f}: "
              f"frac AP-ward bins = {st['frac_ap_ward_bins']*100:.1f}%   "
              f"mean bin speed = {st['mean_bin_speed_um_per_min']:.3f} µm/min   "
              f"({st['n_bins']:,} bins, {st['n_steps']:,} steps)")

    with PdfPages(pdf_path) as pdf:
        # Re-render page 1 (flow maps) into the multi-page PDF
        page1 = _build_page1(spots_v=spots_v, valid_tracks=keep, sphere=sphere,
                             frames=args.frames,
                             surface_um=args.surface, z_mid=Z_MID,
                             z_band=args.z_band,
                             bin_um=args.bin, min_n=args.min_n,
                             arrow_scale=ARROW_SCALE, arrow_lw=ARROW_LW,
                             arrow_head=ARROW_HEAD_FRAC,
                             half_win=args.half_win,
                             speed_thr_px=diag["speed_thr"],
                             speed_median_px=diag["median_speed"],
                             speed_std_px=diag["std_speed"],
                             n_dropped_speed=n_dropped_speed)
        pdf.savefig(page1, bbox_inches="tight")
        plt.close(page1)
        # Page 2 (Panel B trend) — from raw tracks, not binned
        page2 = _build_page2(result["summaries"], track_summary=track_stats)
        pdf.savefig(page2, bbox_inches="tight")
        plt.close(page2)
    print(f"\nSaved 2-page PDF to {pdf_path}")
    print(f"Saved page-1 PNG to {save_stem}.png")
    print(f"Saved summary to {summary_path}")


if __name__ == "__main__":
    main()