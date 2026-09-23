# Vector flow maps — medaka (mk2508) and zebrafish (zb1105)

Binned vector-flow field analyses derived from Stephane's oriented-and-filtered
ultrack tracks for two reference recordings: the medaka `mk2508` low-res
recording (Aug 2025) and the zebrafish `zb1105` recording.

These flow maps were produced in **September 2026** and are deliberately
isolated here from the larger comparative-gastrulation manuscript analysis.

## Contents

```
vector_flow_maps/
├── README.md                 # this file
├── renv.lock                 # R environment lockfile (renv)
├── .gitignore
│
├── scripts/                  # R scripts that produce the flow maps
│   ├── flow_maps_mk2508_backward.R   # mk2508, backward-window variant (5 contraction + 3 internalisation landmarks)
│   ├── flow_maps_mk2508_forward.R    # mk2508, forward-window variant (3 internalisation landmarks)
│   └── flow_maps_zb1105_forward.R    # zb1105, forward-window variant (3 internalisation landmarks)
│
├── inputs/                   # raw oriented+filtered track CSVs (one per dataset)
│   └── tracks_for_flow_maps/
│       ├── tracks_mk_2508_stephane.csv                  # raw mk2508 tracks (Stephane)
│       ├── tracks_mk_2508_stephane_oriented_filtered.csv # mk2508 oriented + filtered (used by mk2508 scripts)
│       ├── tracks_mk_3005_stephane.csv                  # mk3005 tracks — CURRENTLY UNUSED by these scripts
│       ├── tracks_zb_1105.csv                            # raw zb1105 tracks
│       └── tracks_zb_1105_oriented_filtered.csv         # zb1105 oriented + filtered (used by zb1105 script)
│
└── outputs/                  # all per-panel PDF/PNG/CSVs produced by the scripts
    ├── mk2508/
    │   ├── backward_10min/        # 10-min backward-window panels (2-min velocity lag)
    │   ├── backward_30min/        # 30-min backward-window panels (default, 4-frame velocity lag)
    │   ├── forward_internalisation/  # 30-min forward-window panels (internalisation anchors only)
    │   ├── flow_fields.csv        # raw bin-level flow fields (last 30-min backward run)
    │   ├── panel_majority_movement.csv
    │   ├── panel_metadata.csv
    │   └── summary_references.txt  # audit log: voxel size, FI, reference frames, window size
    └── zb1105/
        ├── forward_10min/         # 10-min forward-window panels
        ├── forward_30min/         # 30-min forward-window panels
        ├── 04a_flow_30min_after_internalization_zb1105.pdf  # the 3-panel PDF (root output)
        ├── 04c_majority_movement_zb1105.pdf
        ├── flow_fields.csv
        ├── panel_majority_movement.csv
        ├── panel_metadata.csv
        └── summary_references.txt
```

## How to run

From this directory:

```r
# mk2508, 30-minute backward windows (5 contraction-anchored + 3 internalisation-anchored landmarks)
Rscript scripts/flow_maps_mk2508_backward.R

# mk2508, 30-minute forward windows (3 internalisation-anchored landmarks only)
Rscript scripts/flow_maps_mk2508_forward.R

# zb1105, 30-minute forward windows (3 internalisation-anchored landmarks only)
Rscript scripts/flow_maps_zb1105_forward.R
```

Each script's header comment documents its inputs, outputs, and the
landmarks/voxel size it expects.

## Reference frames (audit log)

### mk2508 (Sep 2026, 0.64 µm isotropic voxels, FI = 30 s)

Anchored to **first contraction** (5 landmarks, absolute frame units):
| t (frame) | landmark                       |
|---:|---|
| 116 | 3 h before first contraction   |
| 160 | start of internalisation       |
| 236 | 2 h before first contraction   |
| 356 | 1 h before first contraction   |
| 476 | first contraction / bulge       |

Anchored to **start of internalisation** (3 landmarks):
| t (frame) | landmark                       |
|---:|---|
|  40 | 1 h before internalisation     |
| 160 | start of internalisation       |
| 280 | 1 h after internalisation      |

### zb1105 (Sep 2026, 1.25 µm isotropic voxels, FI = 120 s)

Anchored to **start of internalisation** (3 landmarks):
| t (frame) | landmark                       |
|---:|---|
|   0 | 1 h before internalisation     |
|  30 | start of internalisation       |
|  60 | 1 h after internalisation      |

## Provenance

These scripts and outputs were originally scattered across
`analysis_output_mk2508/`, `analysis_output_zb1105/`, and
`tracks_for_flow_maps/` at the repository root. The folder structure here
collapses them into one self-contained package so the work can be reviewed
independently of the larger comparative-gastrulation manuscript analysis.
