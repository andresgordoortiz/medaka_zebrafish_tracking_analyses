#!/usr/bin/env python3
"""
Export ultrack segmented nuclei + tracks for viewing in Mastodon.

Converts ultrack's results/ directory into Mastodon-ready files:
  1. segments.zarr  →  BDV HDF5/XML  (so Mastodon displays the 3D nuclei)
  2. tracks.csv     →  (TrackMate XML already exists as ultrack_tracks.xml)

Run this script ON THE CLUSTER where the ultrack results live.

PREREQUISITES:
  pip install pybdv zarr numpy pandas tqdm

USAGE:

  # Full pipeline: convert segments.zarr → BDV + verify tracks XML
  python export_to_mastodon.py /path/to/medaka_25082025

  # Specify voxel size explicitly (Z Y X in µm)
  python export_to_mastodon.py /path/to/medaka_25082025 \
      --voxel-size 1.05 1.05 1.05 --time-interval 90

  # Only convert a subset of timepoints (e.g. first 50)
  python export_to_mastodon.py /path/to/medaka_25082025 \
      --t-start 0 --t-stop 50

  # Skip BDV conversion, just generate TrackMate XML from tracks.csv
  python export_to_mastodon.py /path/to/medaka_25082025 --tracks-only

MASTODON WORKFLOW (on your local machine):
  1. Copy the output/ folder to your local machine
  2. Open Fiji → Plugins → Mastodon
  3. "New Mastodon project" → browse to output/segments_bdv.xml
  4. File → Import → TrackMate file → browse to results/ultrack_tracks.xml
     (or output/tracks_trackmate.xml if regenerated)
  5. You'll see the segmented 3D nuclei with tracking links overlaid.
"""

from __future__ import annotations

import argparse
import sys
import time as time_module
from pathlib import Path
from xml.dom.minidom import parseString
from xml.etree.ElementTree import Element, SubElement, tostring

import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────────────
#  Convert segments.zarr → BDV HDF5
# ──────────────────────────────────────────────────────────────────────

def convert_segments_to_bdv(
    segments_zarr_path: Path,
    out_path: Path,
    voxel_size: tuple[float, float, float],
    time_interval: float,
    t_start: int | None = None,
    t_stop: int | None = None,
    downscale_factors: list[list[int]] | None = None,
    n_workers: int | None = None,
) -> Path:
    """
    Read ultrack's segments.zarr (4D label array, TZYX) and write it as
    a BDV HDF5 + XML file that Mastodon can open.

    Each voxel in segments.zarr contains the track_id of the nucleus
    occupying that voxel (0 = background).  Mastodon will render these
    as an intensity image — nuclei appear bright on a dark background.

    Parameters
    ----------
    segments_zarr_path : Path to results/segments.zarr
    out_path : Path for output .h5 file (the .xml is produced alongside)
    voxel_size : (dz, dy, dx) in µm
    time_interval : seconds between frames
    t_start, t_stop : optional timepoint range (inclusive start, exclusive stop)
    downscale_factors : multi-resolution pyramid factors
    n_workers : number of parallel workers (default: all CPUs)

    Returns
    -------
    Path to the .xml file
    """
    try:
        import pybdv  # noqa: F401 — verify it's available before starting
    except ImportError:
        print("ERROR: pybdv is required.  Install with:  pip install pybdv")
        sys.exit(1)

    import zarr

    print(f"Opening {segments_zarr_path} ...")
    store = zarr.open(str(segments_zarr_path), mode="r")

    # Detect structure: single 4D array vs group of 3D arrays
    if isinstance(store, zarr.Array):
        # Single 4D array (T, Z, Y, X)
        n_timepoints = store.shape[0]
        spatial_shape = store.shape[1:]
        is_4d = True
        print(f"  4D array: shape={store.shape}, dtype={store.dtype}")
    elif isinstance(store, zarr.Group):
        # Group with per-timepoint 3D arrays keyed by integer
        keys = sorted([int(k) for k in store.keys()])
        n_timepoints = len(keys)
        first = np.asarray(store[str(keys[0])])
        spatial_shape = first.shape
        is_4d = False
        print(f"  Zarr group: {n_timepoints} timepoints, spatial shape={spatial_shape}")
    else:
        raise ValueError(f"Unexpected zarr structure at {segments_zarr_path}")

    if t_start is None:
        t_start = 0
    if t_stop is None:
        t_stop = n_timepoints
    t_stop = min(t_stop, n_timepoints)
    print(f"  Converting timepoints {t_start}–{t_stop - 1} "
          f"({t_stop - t_start} frames)")
    print(f"  Voxel size: {voxel_size} µm,  time interval: {time_interval} s")

    if downscale_factors is None:
        # Reasonable pyramid for typical SPIM data
        downscale_factors = [[2, 2, 2], [2, 2, 2], [2, 2, 2]]

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Remove previous output to avoid "already present" warnings
    if out_path.exists():
        out_path.unlink()
    xml_candidate = out_path.with_suffix(".xml")
    if xml_candidate.exists():
        xml_candidate.unlink()

    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed
    if n_workers is None:
        n_workers = os.cpu_count() or 1

    timepoints = list(range(t_start, t_stop))
    n_tp = len(timepoints)

    # Write each timepoint to its own temp HDF5 in parallel,
    # then merge into the final file.
    tmp_dir = out_path.parent / "_bdv_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Using {n_workers} parallel workers")

    wall_t0 = time_module.time()

    # Submit parallel jobs
    futures = {}
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for t in timepoints:
            tmp_h5 = tmp_dir / f"tp_{t:06d}.h5"
            fut = pool.submit(
                _process_timepoint,
                str(segments_zarr_path),
                str(tmp_h5),
                t,
                0,  # write as timepoint 0 in temp file
                downscale_factors,
                list(voxel_size),
                is_4d,
            )
            futures[fut] = t

        done = 0
        for fut in as_completed(futures):
            done += 1
            t = futures[fut]
            exc = fut.exception()
            if exc:
                print(f"\n  ERROR at timepoint {t}: {exc}")
                raise exc
            if done % 10 == 0 or done == n_tp:
                elapsed = time_module.time() - wall_t0
                rate = elapsed / done
                eta = rate * (n_tp - done)
                print(f"  [{done}/{n_tp}] timepoints processed  "
                      f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

    # Merge temp HDF5 files into the final BDV file
    print("  Merging into final BDV HDF5 ...")
    import h5py
    import shutil
    merge_t0 = time_module.time()

    # The first temp file was written by make_bdv, which also created
    # a matching .xml — copy it as the starting point for the final file.
    first_tmp = tmp_dir / f"tp_{timepoints[0]:06d}.h5"
    first_xml = first_tmp.with_suffix(".xml")
    shutil.copy2(str(first_tmp), str(out_path))
    if first_xml.exists():
        shutil.copy2(str(first_xml), str(out_path.with_suffix(".xml")))

    # Append remaining timepoints into the final HDF5
    # pybdv stores data in nested groups: t00000/s00/0/cells, t00000/s00/1/cells, ...
    # We need to recursively copy the entire group tree.
    def _copy_group(src_group, dst_group):
        """Recursively copy all datasets and subgroups from src to dst."""
        for key in src_group:
            item = src_group[key]
            if isinstance(item, h5py.Dataset):
                ds = dst_group.create_dataset(
                    key, data=item[:],
                    chunks=item.chunks,
                    compression=item.compression,
                    compression_opts=item.compression_opts,
                )
                for k, v in item.attrs.items():
                    ds.attrs[k] = v
            elif isinstance(item, h5py.Group):
                sub = dst_group.create_group(key)
                for k, v in item.attrs.items():
                    sub.attrs[k] = v
                _copy_group(item, sub)

    if len(timepoints) > 1:
        with h5py.File(str(out_path), "a") as dst:
            for i, t in enumerate(timepoints[1:], start=1):
                tp_key = f"t{i:05d}"
                tmp_h5 = tmp_dir / f"tp_{t:06d}.h5"
                with h5py.File(str(tmp_h5), "r") as src:
                    src_tp_grp = src["t00000"]
                    dst_tp_grp = dst.create_group(tp_key)
                    for k, v in src_tp_grp.attrs.items():
                        dst_tp_grp.attrs[k] = v
                    _copy_group(src_tp_grp, dst_tp_grp)

    # Fix XML to declare all timepoints (make_bdv only wrote tp 0)
    _fix_bdv_xml_timepoints(out_path.with_suffix(".xml"), n_tp)

    # Cleanup temp files
    shutil.rmtree(str(tmp_dir), ignore_errors=True)

    elapsed = time_module.time() - wall_t0
    merge_elapsed = time_module.time() - merge_t0
    xml_path = out_path.with_suffix(".xml")
    print(f"\nBDV conversion done in {elapsed:.0f}s "
          f"(merge: {merge_elapsed:.0f}s)")
    print(f"  HDF5: {out_path}")
    print(f"  XML:  {xml_path}")
    return xml_path


def _remap_to_uint16(vol: np.ndarray) -> np.ndarray:
    """Remap label volume to contiguous 1..N uint16 values."""
    unique_ids = np.unique(vol)
    unique_ids = unique_ids[unique_ids != 0]
    if len(unique_ids) > 0:
        lut = np.zeros(int(vol.max()) + 1, dtype=np.uint16)
        lut[unique_ids] = np.arange(1, len(unique_ids) + 1, dtype=np.uint16)
        return lut[vol]
    return vol.astype(np.uint16)


def _process_timepoint(
    zarr_path: str,
    out_h5: str,
    t: int,
    tp_index: int,
    downscale_factors: list,
    voxel_size: list,
    is_4d: bool,
) -> None:
    """Read one timepoint from zarr, remap labels, write to temp HDF5."""
    import zarr as zarr_mod
    from pybdv import make_bdv

    store = zarr_mod.open(zarr_path, mode="r")
    if is_4d:
        vol = np.asarray(store[t])
    else:
        vol = np.asarray(store[str(t)])

    vol = _remap_to_uint16(vol)

    make_bdv(
        vol,
        out_h5,
        setup_id=0,
        timepoint=tp_index,
        downscale_factors=downscale_factors,
        downscale_mode="nearest",
        resolution=voxel_size,
        unit="micrometer",
    )


def _fix_bdv_xml_timepoints(xml_path: Path, n_timepoints: int) -> None:
    """Update BDV XML to declare all timepoints and view registrations."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(str(xml_path))
    root = tree.getroot()

    # Find or create SequenceDescription/Timepoints
    seq = root.find(".//SequenceDescription")
    if seq is None:
        return
    tp_el = seq.find("Timepoints")
    if tp_el is not None:
        seq.remove(tp_el)

    tp_el = ET.SubElement(seq, "Timepoints", type="range")
    first = ET.SubElement(tp_el, "first")
    first.text = "0"
    last = ET.SubElement(tp_el, "last")
    last.text = str(n_timepoints - 1)

    # Fix HDF5 path: the XML was copied from a temp file so it still
    # references e.g. "tp_000000.h5" — update to the actual .h5 filename.
    h5_name = xml_path.with_suffix(".h5").name
    for hdf5_el in root.iter("hdf5"):
        hdf5_el.text = h5_name
    for il in root.iter("ImageLoader"):
        hdf5_child = il.find("hdf5")
        if hdf5_child is not None:
            hdf5_child.text = h5_name

    # Fix ViewRegistrations: pybdv only wrote registration for tp 0.
    # BDV/Mastodon needs a ViewRegistration for every (timepoint, setup).
    vr_el = root.find(".//ViewRegistrations")
    if vr_el is not None:
        # Collect existing registrations and their transforms
        existing = {}  # (tp, setup) -> list of ViewTransform elements
        for vr in vr_el.findall("ViewRegistration"):
            tp = int(vr.get("timepoint", "0"))
            setup = int(vr.get("setup", "0"))
            existing[(tp, setup)] = list(vr)

        # Determine which setups exist (usually just setup 0)
        setups = sorted(set(s for _, s in existing.keys()))
        if not setups:
            setups = [0]

        # For each setup, get the template from tp 0 and replicate
        for setup in setups:
            template_transforms = existing.get((0, setup), [])
            for tp in range(n_timepoints):
                if (tp, setup) in existing:
                    continue
                vr_new = ET.SubElement(vr_el, "ViewRegistration",
                                       timepoint=str(tp), setup=str(setup))
                for vt in template_transforms:
                    # Deep copy the ViewTransform element
                    vt_new = ET.SubElement(vr_new, vt.tag, **vt.attrib)
                    vt_new.text = vt.text
                    vt_new.tail = vt.tail
                    for child in vt:
                        child_new = ET.SubElement(vt_new, child.tag, **child.attrib)
                        child_new.text = child.text
                        child_new.tail = child.tail

    tree.write(str(xml_path), xml_declaration=True, encoding="utf-8")


# ──────────────────────────────────────────────────────────────────────
#  Read ultrack tracks CSV
# ──────────────────────────────────────────────────────────────────────

def load_ultrack_tracks(path: Path) -> pd.DataFrame:
    """Load ultrack tracks.csv (raw format: track_id, t, z, y, x, parent_track_id)."""
    df = pd.read_csv(path, low_memory=False)

    # Normalise column names — accept both raw ultrack and oriented formats
    col_map = {}
    if "track_id" in df.columns and "TRACK_ID" not in df.columns:
        col_map["track_id"] = "TRACK_ID"
    if "t" in df.columns and "FRAME" not in df.columns:
        col_map["t"] = "FRAME"
    if "z" in df.columns and "POSITION_Z" not in df.columns:
        col_map["z"] = "POSITION_Z"
    if "y" in df.columns and "POSITION_Y" not in df.columns:
        col_map["y"] = "POSITION_Y"
    if "x" in df.columns and "POSITION_X" not in df.columns:
        col_map["x"] = "POSITION_X"
    if col_map:
        df = df.rename(columns=col_map)

    required = ["TRACK_ID", "FRAME", "POSITION_X", "POSITION_Y", "POSITION_Z"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"CSV is missing columns: {missing}")

    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=required)
    df["FRAME"] = df["FRAME"].astype(int)
    df["TRACK_ID"] = df["TRACK_ID"].astype(int)

    # Ensure every row has a unique spot id
    if "id" not in df.columns:
        df["id"] = np.arange(1, len(df) + 1)
    else:
        df["id"] = pd.to_numeric(df["id"], errors="coerce").astype(int)

    # Build parent_id from track continuity if not present
    if "parent_id" not in df.columns:
        df = df.sort_values(["TRACK_ID", "FRAME"]).reset_index(drop=True)
        df["parent_id"] = -1
        for _, grp in df.groupby("TRACK_ID"):
            idxs = grp.index.tolist()
            for i in range(1, len(idxs)):
                df.loc[idxs[i], "parent_id"] = int(df.loc[idxs[i - 1], "id"])
    else:
        df["parent_id"] = (pd.to_numeric(df["parent_id"], errors="coerce")
                           .fillna(-1).astype(int))

    return df


# ──────────────────────────────────────────────────────────────────────
#  Export: TrackMate XML  (Mastodon can import directly)
# ──────────────────────────────────────────────────────────────────────

def export_trackmate_xml(
    df: pd.DataFrame,
    out_path: Path,
    radius: float,
    pixel_width: float,
    pixel_height: float,
    voxel_depth: float,
    time_interval: float,
) -> None:
    """
    Write a TrackMate-compatible XML from a tracks DataFrame.

    Reproduces ultrack's tracks_layer_to_trackmate() without needing
    ultrack installed.
    """
    df = df.sort_values(["TRACK_ID", "FRAME"]).reset_index(drop=True)
    n_frames = int(df["FRAME"].max()) + 1

    root = Element("TrackMate", version="7.11.1")
    model = SubElement(root, "Model", spatialunits="pixels", timeunits="frames")

    # Feature declarations
    feat_decl = SubElement(model, "FeatureDeclarations")
    spot_feats = SubElement(feat_decl, "SpotFeatures")
    for fname, sname, dim in [
        ("QUALITY", "Quality", "QUALITY"),
        ("POSITION_X", "X", "POSITION"),
        ("POSITION_Y", "Y", "POSITION"),
        ("POSITION_Z", "Z", "POSITION"),
        ("FRAME", "Frame", "NONE"),
        ("RADIUS", "Radius", "LENGTH"),
    ]:
        SubElement(spot_feats, "Feature", feature=fname, name=sname,
                   shortname=sname, dimension=dim, isint="false")
    SubElement(feat_decl, "EdgeFeatures")
    SubElement(feat_decl, "TrackFeatures")

    # AllSpots
    all_spots = SubElement(model, "AllSpots", nspots=str(len(df)))
    for frame, grp in df.groupby("FRAME"):
        sif = SubElement(all_spots, "SpotsInFrame", frame=str(int(frame)))
        for _, row in grp.iterrows():
            SubElement(sif, "Spot",
                       ID=str(int(row["id"])),
                       name=str(int(row["id"])),
                       QUALITY="1.0",
                       VISIBILITY="1",
                       FRAME=str(int(row["FRAME"])),
                       RADIUS=f"{radius:.1f}",
                       POSITION_X=f"{row['POSITION_X']:.4f}",
                       POSITION_Y=f"{row['POSITION_Y']:.4f}",
                       POSITION_Z=f"{row['POSITION_Z']:.4f}")

    # AllTracks
    all_tracks = SubElement(model, "AllTracks")
    filtered_tracks = SubElement(model, "FilteredTracks")

    for tid, grp in df.groupby("TRACK_ID"):
        grp = grp.sort_values("FRAME")
        edges = []
        for i in range(1, len(grp)):
            src_id = int(grp.iloc[i - 1]["id"])
            tgt_id = int(grp.iloc[i]["id"])
            edge_t = (grp.iloc[i - 1]["FRAME"] + grp.iloc[i]["FRAME"]) / 2.0
            edges.append((src_id, tgt_id, edge_t))

        # Handle division links via parent_track_id
        first_row = grp.iloc[0]
        if int(first_row["parent_id"]) > 0:
            src_id = int(first_row["parent_id"])
            tgt_id = int(first_row["id"])
            edge_t = float(first_row["FRAME"]) - 0.5
            edges.insert(0, (src_id, tgt_id, edge_t))

        if not edges:
            continue

        track_el = SubElement(all_tracks, "Track",
                              TRACK_ID=str(int(tid)),
                              NUMBER_SPOTS=str(len(grp)),
                              NUMBER_GAPS="0",
                              TRACK_START=str(int(grp["FRAME"].min())),
                              TRACK_STOP=str(int(grp["FRAME"].max())),
                              name=f"Track_{tid}")
        for src_id, tgt_id, edge_t in edges:
            SubElement(track_el, "Edge",
                       SPOT_SOURCE_ID=str(src_id),
                       SPOT_TARGET_ID=str(tgt_id),
                       EDGE_TIME=f"{edge_t:.1f}")

        SubElement(filtered_tracks, "TrackID", TRACK_ID=str(int(tid)))

    # Settings (image metadata)
    settings = SubElement(root, "Settings")
    SubElement(settings, "InitialSpotFilter", feature="QUALITY",
               value="0.0", isabove="true")
    SubElement(settings, "SpotFilterCollection")
    SubElement(settings, "TrackFilterCollection")
    SubElement(settings, "ImageData",
               filename="None", folder="None",
               width="0", height="0", depth="0",
               nslices="1", nframes=str(n_frames),
               pixelwidth=str(pixel_width),
               pixelheight=str(pixel_height),
               voxeldepth=str(voxel_depth),
               timeinterval=str(time_interval))

    xml_str = parseString(tostring(root, encoding="unicode")).toprettyxml(indent="  ")
    out_path.write_text(xml_str, encoding="utf-8")

    print(f"TrackMate XML written → {out_path}")
    print(f"  {len(df)} spots, {df['TRACK_ID'].nunique()} tracks, {n_frames} frames")


# ──────────────────────────────────────────────────────────────────────
#  CLI
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Export ultrack results (segmented nuclei + tracks) "
                    "for viewing in Mastodon.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "ultrack_dir", type=Path,
        help="Root ultrack directory containing results/ "
             "(e.g. /path/to/medaka_25082025)")
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=None,
        help="Output directory (default: <ultrack_dir>/mastodon_export)")
    parser.add_argument(
        "--voxel-size", type=float, nargs=3, default=[1.0, 1.0, 1.0],
        metavar=("Z", "Y", "X"),
        help="Voxel size in µm (Z Y X).  Default: 1.0 1.0 1.0")
    parser.add_argument(
        "--time-interval", type=float, default=30.0,
        help="Time interval between frames in seconds (default: 30.0)")
    parser.add_argument(
        "--radius", type=float, default=5.0,
        help="Default spot radius in pixels for TrackMate XML (default: 5.0)")
    parser.add_argument(
        "--t-start", type=int, default=None,
        help="First timepoint to convert (default: 0)")
    parser.add_argument(
        "--t-stop", type=int, default=None,
        help="Last timepoint (exclusive) to convert (default: all)")
    parser.add_argument(
        "--tracks-only", action="store_true",
        help="Skip BDV conversion; only regenerate TrackMate XML from tracks.csv")
    parser.add_argument(
        "--downscale", type=str, default=None,
        help="Downscale factors as JSON, e.g. '[[2,2,2],[2,2,2],[2,2,2]]'")
    parser.add_argument(
        "--workers", type=int, default=None,
        help="Number of parallel workers for BDV conversion (default: all CPUs)")

    args = parser.parse_args()

    ultrack_dir = args.ultrack_dir.resolve()
    results_dir = ultrack_dir / "results"

    # Validate directory structure
    if not results_dir.is_dir():
        print(f"ERROR: {results_dir} not found.", file=sys.stderr)
        print(f"Expected directory layout:", file=sys.stderr)
        print(f"  {ultrack_dir}/", file=sys.stderr)
        print(f"    results/", file=sys.stderr)
        print(f"      segments.zarr/", file=sys.stderr)
        print(f"      tracks.csv", file=sys.stderr)
        print(f"      ultrack_tracks.xml  (optional)", file=sys.stderr)
        sys.exit(1)

    segments_path = results_dir / "segments.zarr"
    tracks_csv_path = results_dir / "tracks.csv"
    existing_xml = results_dir / "ultrack_tracks.xml"

    out_dir = args.output_dir or ultrack_dir / "mastodon_export"
    out_dir.mkdir(parents=True, exist_ok=True)

    voxel_size = tuple(args.voxel_size)
    downscale_factors = None
    if args.downscale:
        import json
        downscale_factors = json.loads(args.downscale)

    print("=" * 60)
    print("  ultrack → Mastodon exporter")
    print("=" * 60)
    print(f"  ultrack dir:    {ultrack_dir}")
    print(f"  output dir:     {out_dir}")
    print(f"  voxel size:     {voxel_size} µm")
    print(f"  time interval:  {args.time_interval} s")
    print()

    # ── Step 1: Convert segments.zarr → BDV HDF5 ──
    bdv_xml = None
    if not args.tracks_only:
        if not segments_path.is_dir():
            print(f"WARNING: {segments_path} not found — skipping BDV conversion.")
            print("         Use --tracks-only to only export tracks.")
        else:
            bdv_h5 = out_dir / "segments_bdv.h5"
            bdv_xml = convert_segments_to_bdv(
                segments_path,
                bdv_h5,
                voxel_size=voxel_size,
                time_interval=args.time_interval,
                t_start=args.t_start,
                t_stop=args.t_stop,
                downscale_factors=downscale_factors,
                n_workers=args.workers,
            )
            print()

    # ── Step 2: Tracks ──
    # Check if ultrack already produced a TrackMate XML
    if existing_xml.exists():
        print(f"TrackMate XML already exists: {existing_xml}")
        print("  You can import this directly in Mastodon.")
    elif tracks_csv_path.exists():
        print(f"No ultrack_tracks.xml found, generating from {tracks_csv_path} ...")
        df = load_ultrack_tracks(tracks_csv_path)
        xml_out = out_dir / "tracks_trackmate.xml"
        export_trackmate_xml(
            df, xml_out,
            radius=args.radius,
            pixel_width=voxel_size[2],
            pixel_height=voxel_size[1],
            voxel_depth=voxel_size[0],
            time_interval=args.time_interval,
        )
    else:
        print(f"WARNING: No tracks found at {tracks_csv_path} or {existing_xml}")

    # Also regenerate from CSV if explicitly requested
    if args.tracks_only and tracks_csv_path.exists():
        print(f"\nRegenerating TrackMate XML from {tracks_csv_path} ...")
        df = load_ultrack_tracks(tracks_csv_path)
        xml_out = out_dir / "tracks_trackmate.xml"
        export_trackmate_xml(
            df, xml_out,
            radius=args.radius,
            pixel_width=voxel_size[2],
            pixel_height=voxel_size[1],
            voxel_depth=voxel_size[0],
            time_interval=args.time_interval,
        )

    # ── Summary ──
    print()
    print("=" * 60)
    print("  DONE — How to view in Mastodon:")
    print("=" * 60)
    print()
    if bdv_xml:
        print(f"  1. Copy {out_dir}/ to your local machine")
        print(f"  2. Open Fiji → Plugins → Mastodon")
        print(f"  3. 'New Mastodon project' → {bdv_xml.name}")
        print(f"  4. File → Import → TrackMate file →", end=" ")
        if existing_xml.exists():
            print(existing_xml.name)
        else:
            print("tracks_trackmate.xml")
        print()
        print("  You'll see the segmented 3D nuclei with tracking links.")
    else:
        print("  BDV image was not generated (--tracks-only or segments.zarr missing).")
        print("  You can still import the TrackMate XML into an existing Mastodon project.")


if __name__ == "__main__":
    main()
