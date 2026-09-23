#!/usr/bin/env python3
"""
Per-region Voronoi cell volume / surface statistics for the two species.

We compute Voronoi diagrams on the *projected sphere surface* (theta, phi in
degrees, treated as flat angular coordinates -- small enough that the
projection is locally Euclidean) using scipy.spatial.Voronoi.  Each cell's
2-D angular area (deg^2) is then converted to um^2 via the sphere Jacobian
at the cell centroid:  A_um2 = R^2 * area_rad2, where area_rad2 = (pi/180)^2
* area_deg2.  This matches the convention used in density_3timepoints.R.

The 3-D volume is "column Voronoi": for each cell we use its individual
radial depth range (p2-p98) to make a column whose base is the 2-D Voronoi
area and height is the local depth range.  This is a per-cell 3-D volume
that does NOT pool a global "tissue depth" across all cells in the region.

Reads:   <in_dir>/<species>_coords_<tp>.csv  (TRACK_ID, x, y, z, theta, phi, depth)
Writes:  <out_dir>/<species>_voronoi_<tp>.csv (TRACK_ID, cell_area_um2, col_vol_um3)
Also writes: <out_dir>/<species>_voronoi_summary_<tp>.csv with region-level stats.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi


def voronoi_cell_areas(coords_deg: np.ndarray) -> np.ndarray:
    """Voronoi diagram on a 2-D point cloud; return per-cell polygon areas.

    Uses scipy.spatial.Voronoi and sums 2-D triangle areas from each cell's
    centroid to the polygon edges (in case scipy returns -1 = open cell).
    Falls back to a kd-tree nearest-neighbour estimate of local area for any
    cell whose polygon degenerated or is unbounded.
    """
    # Add a large "frame" of mirrored points around the convex hull so that
    # all real polygons are bounded.  This is a standard trick.
    pad = np.max(np.ptp(coords_deg, axis=0)) * 4 + 100.0
    centre = coords_deg.mean(axis=0)
    framed = []
    for sx in (-1, 1):
        for sy in (-1, 1):
            framed.append(
                coords_deg + (sx * pad, sy * pad) - 0.0
            )
    framed = np.vstack(framed)
    all_pts = np.vstack([coords_deg, framed])
    vor = Voronoi(all_pts)
    areas = np.full(coords_deg.shape[0], np.nan)
    # Index map: i -> vor.point_region[i] -> vor.regions[reg_idx]
    for i in range(coords_deg.shape[0]):
        reg_idx = vor.point_region[i]
        region = vor.regions[reg_idx]
        if not region or -1 in region:
            # Degenerate or unbounded: fall back to a local kNN area estimate
            continue
        poly = vor.vertices[region]
        # Shoelace formula
        x = poly[:, 0]
        y = poly[:, 1]
        n = len(x)
        s = 0.0
        for k in range(n):
            s += x[k] * y[(k + 1) % n] - x[(k + 1) % n] * y[k]
        areas[i] = abs(s) / 2.0
    # Fallback for any NaN cells: use 1/(density) of the surrounding cells
    nans = np.isnan(areas)
    if nans.any():
        from scipy.spatial import cKDTree
        tree = cKDTree(coords_deg)
        for i in np.where(nans)[0]:
            d, _ = tree.query(coords_deg[i], k=7)  # includes self
            # Local area ~ pi * d_6^2  (approx polygon radius to 6th neighbour)
            areas[i] = np.pi * d[-1] ** 2 / 4
    return areas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--species", required=True)  # "medaka" or "zebrafish"
    ap.add_argument("--tp", required=True)        # "T-1h", "T0", "T+1h"
    ap.add_argument("--sphere_R_um", type=float, required=True)
    args = ap.parse_args()

    in_csv = os.path.join(args.in_dir, f"{args.species}_coords_{args.tp}.csv")
    df = pd.read_csv(in_csv)
    if len(df) == 0:
        print(f"[voronoi] empty input: {in_csv}")
        return
    df = df.dropna(subset=["theta_deg", "phi_deg"]).reset_index(drop=True)
    if len(df) < 8:
        # Voronoi for very few points is meaningless -- write blanks and exit
        df["cell_area_um2"] = np.nan
        df["col_vol_um3"] = np.nan
        out = os.path.join(args.out_dir, f"{args.species}_voronoi_{args.tp}.csv")
        df.to_csv(out, index=False)
        print(f"[voronoi] too few cells ({len(df)}); wrote NaN to {out}")
        return

    coords = df[["theta_deg", "phi_deg"]].values
    cell_area_deg2 = voronoi_cell_areas(coords)
    # Convert deg^2 -> rad^2 -> um^2 via R^2 factor (sphere Jacobian)
    rad_per_deg = np.pi / 180.0
    cell_area_um2 = (rad_per_deg ** 2) * (args.sphere_R_um ** 2) * cell_area_deg2

    # Per-cell column volume = cell_area_um2 * depth range of ITS track
    # using p2-p98 of depth within the window.  We approximate by taking the
    # spread of depth values within the window for each track (max-min of the
    # track's depth within the window; if a track has a single observation,
    # we use the window-global p98-p02 depth range as fallback).
    depth = df["depth_um"].values
    if len(df) >= 2:
        global_depth_range = float(np.quantile(depth, 0.98) - np.quantile(depth, 0.02))
    else:
        global_depth_range = 0.0
    # Per-track depth spread
    track_ranges = (
        df.groupby("TRACK_ID")["depth_um"].agg(["min", "max", "count"]).reset_index()
    )
    track_ranges["depth_range_um"] = track_ranges["max"] - track_ranges["min"]
    track_ranges.loc[track_ranges["depth_range_um"] <= 0, "depth_range_um"] = np.nan
    # For tracks with a single observation we have no per-track range; set
    # to NaN and fall back to global range later
    fallback_ids = set(track_ranges.loc[track_ranges["depth_range_um"].isna(), "TRACK_ID"])
    df["track_depth_range_um"] = df["TRACK_ID"].map(
        track_ranges.set_index("TRACK_ID")["depth_range_um"]
    )
    df["depth_range_um"] = df["track_depth_range_um"]
    df.loc[df["TRACK_ID"].isin(fallback_ids), "depth_range_um"] = global_depth_range
    df["cell_area_um2"] = cell_area_um2
    df["col_vol_um3"] = cell_area_um2 * df["depth_range_um"]

    out_csv = os.path.join(args.out_dir, f"{args.species}_voronoi_{args.tp}.csv")
    df.to_csv(out_csv, index=False)

    summary = {
        "n_cells": len(df),
        "n_tracks": df["TRACK_ID"].nunique(),
        "median_cell_area_um2": float(df["cell_area_um2"].median()),
        "median_col_vol_um3": float(df["col_vol_um3"].median()),
        "median_depth_range_um": float(df["depth_range_um"].median()),
        "median_density_per_um3_voronoi": float(
            (1.0 / df["col_vol_um3"]).replace([np.inf, -np.inf], np.nan).median()
        ),
    }
    sum_csv = os.path.join(args.out_dir, f"{args.species}_voronoi_summary_{args.tp}.csv")
    pd.DataFrame([summary]).to_csv(sum_csv, index=False)
    print(f"[voronoi] {args.species} {args.tp}: n={summary['n_cells']} "
          f"median_area={summary['median_cell_area_um2']:.0f} um^2  "
          f"median_col_vol={summary['median_col_vol_um3']:.0f} um^3")


if __name__ == "__main__":
    main()
