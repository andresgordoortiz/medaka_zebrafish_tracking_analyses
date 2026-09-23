#!/usr/bin/env python3
"""
3D nuclei + tracks viewer — see actual segmented nuclei across time in napari.

Loads ultrack's segments.zarr as a Labels layer (lazy via dask) so you see
the real 3D shapes of every nucleus, plus tracks.csv as a Tracks layer for
lineage / motion visualisation.

USAGE:
  # Minimal: just the segmentation
  python nuclei_tracks_viewer.py /path/to/results/segments.zarr

  # With tracks overlay
  python nuclei_tracks_viewer.py /path/to/results/segments.zarr \
      --tracks /path/to/results/tracks.csv

  # Specify voxel size (Z Y X in µm) and time interval
  python nuclei_tracks_viewer.py /path/to/results/segments.zarr \
      --tracks /path/to/results/tracks.csv \
      --voxel-size 1.05 1.05 1.05 --time-interval 30
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_segments_lazy(zarr_path: Path) -> tuple:
    """
    Open segments.zarr and return a dask array (TZYX) + metadata.

    Supports both layouts:
      - Single 4D zarr array  (shape = T, Z, Y, X)
      - Zarr group with per-timepoint 3D arrays keyed by "0", "1", ...

    Works with zarr v2 and v3.
    """
    import dask.array as da
    import zarr

    store = zarr.open(str(zarr_path), mode="r")

    if isinstance(store, zarr.Array):
        # Already a 4D array — wrap in dask for lazy loading (path-based
        # to stay compatible with zarr v2 and v3)
        data = da.from_zarr(str(zarr_path))
        n_tp = data.shape[0]
        spatial_shape = data.shape[1:]
        print(f"  4D zarr array: shape={data.shape}, dtype={data.dtype}")
        return data, n_tp, spatial_shape

    if isinstance(store, zarr.Group):
        keys = sorted(int(k) for k in store.keys())
        n_tp = len(keys)
        first = store[str(keys[0])]
        spatial_shape = first.shape
        dtype = first.dtype

        # Stack per-timepoint arrays lazily (use path + component for compat)
        lazy_frames = []
        for k in keys:
            lazy_frames.append(
                da.from_zarr(str(zarr_path), component=str(k))
            )
        data = da.stack(lazy_frames, axis=0)
        print(f"  Zarr group: {n_tp} timepoints, spatial={spatial_shape}, "
              f"dtype={dtype}")
        return data, n_tp, spatial_shape

    raise ValueError(f"Unexpected zarr layout at {zarr_path}")


def load_tracks(tracks_path: Path) -> np.ndarray:
    """
    Load ultrack tracks.csv and return an array shaped for napari Tracks layer.

    napari Tracks format: (N, 4) array with columns [track_id, t, z, y, x].
    """
    df = pd.read_csv(tracks_path, low_memory=False)

    # Normalise column names
    col_map = {}
    if "track_id" in df.columns:
        col_map["track_id"] = "TRACK_ID"
    if "t" in df.columns:
        col_map["t"] = "FRAME"
    if "z" in df.columns:
        col_map["z"] = "POSITION_Z"
    if "y" in df.columns:
        col_map["y"] = "POSITION_Y"
    if "x" in df.columns:
        col_map["x"] = "POSITION_X"
    if col_map:
        df = df.rename(columns=col_map)

    required = ["TRACK_ID", "FRAME", "POSITION_Z", "POSITION_Y", "POSITION_X"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"tracks CSV missing columns: {missing}")

    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=required)
    df = df.sort_values(["TRACK_ID", "FRAME"]).reset_index(drop=True)

    tracks_arr = df[["TRACK_ID", "FRAME",
                     "POSITION_Z", "POSITION_Y", "POSITION_X"]].to_numpy()

    n_tracks = df["TRACK_ID"].nunique()
    n_spots = len(df)
    n_frames = int(df["FRAME"].max()) - int(df["FRAME"].min()) + 1
    print(f"  {n_spots} spots, {n_tracks} tracks, {n_frames} frames")

    return tracks_arr


def main():
    parser = argparse.ArgumentParser(
        description="View 3D segmented nuclei + tracks in napari.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "segments_zarr", type=Path,
        help="Path to segments.zarr (ultrack output, 4D label volume)")
    parser.add_argument(
        "--tracks", type=Path, default=None,
        help="Path to tracks.csv (ultrack format: track_id, t, z, y, x)")
    parser.add_argument(
        "--voxel-size", type=float, nargs=3, default=[1.0, 1.0, 1.0],
        metavar=("Z", "Y", "X"),
        help="Voxel size in µm (Z Y X). Default: 1.0 1.0 1.0")
    parser.add_argument(
        "--time-interval", type=float, default=30.0,
        help="Time between frames in seconds (default: 30)")
    parser.add_argument(
        "--3d", dest="start_3d", action="store_true",
        help="Start in 3D rendering mode (default: slice view)")

    args = parser.parse_args()

    # ── Wayland / DPI fixes ──
    if os.environ.get("XDG_SESSION_TYPE") == "wayland":
        os.environ["QT_QPA_PLATFORM"] = "wayland"
    os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")

    import napari

    zarr_path = args.segments_zarr.resolve()
    voxel_size = tuple(args.voxel_size)  # (Z, Y, X)
    scale_4d = (1,) + voxel_size         # (T, Z, Y, X)

    # ── Load segmentation ──
    print(f"Loading segments: {zarr_path}")
    seg_data, n_tp, spatial_shape = load_segments_lazy(zarr_path)

    # ── Create napari viewer ──
    viewer = napari.Viewer(
        title=f"Nuclei viewer — {zarr_path.parent.name}",
        ndisplay=3 if args.start_3d else 2,
    )

    # Labels layer — actual 3D nuclei at each timepoint
    labels_layer = viewer.add_labels(
        seg_data,
        name="nuclei",
        scale=scale_4d,
        opacity=0.7,
    )

    print(f"  Labels layer added: {n_tp} timepoints, "
          f"spatial {spatial_shape}, scale={scale_4d}")

    # ── Load tracks ──
    if args.tracks is not None:
        tracks_path = args.tracks.resolve()
        print(f"Loading tracks: {tracks_path}")
        tracks_arr = load_tracks(tracks_path)

        # Scale track coordinates to physical units
        # tracks_arr columns: [track_id, t, z, y, x]
        tracks_scaled = tracks_arr.copy()
        tracks_scaled[:, 2] *= voxel_size[0]  # Z
        tracks_scaled[:, 3] *= voxel_size[1]  # Y
        tracks_scaled[:, 4] *= voxel_size[2]  # X

        tracks_layer = viewer.add_tracks(
            tracks_scaled,
            name="tracks",
            tail_width=2,
            tail_length=30,
            head_length=0,
            color_by="track_id",
        )
        print("  Tracks layer added")

    # ── Configure viewer ──
    # Set the time axis as the slider dimension
    viewer.dims.axis_labels = ["t", "z", "y", "x"]

    # Set a reasonable initial timepoint
    viewer.dims.set_point(0, 0)

    print(f"\nReady — {n_tp} timepoints, scrub the 't' slider to browse.")
    print("  Press '3' to toggle 3D / slice view")
    print("  Ctrl+Shift+E to toggle labels visibility")

    napari.run()


if __name__ == "__main__":
    main()
