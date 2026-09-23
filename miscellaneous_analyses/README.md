# Comparative gastrulation manuscript — miscellaneous analyses

This package contains the comparative **medaka vs. zebrafish** gastrulation
manuscript analysis: 4D lightsheet whole-embryo nuclei tracking, per-timepoint
nuclear statistics, track-level gastrulation dynamics, density and margin
analyses, and a oneshot comparison across the two species.

The recent (Sep 2026) **vector flow maps** work for `mk2508` and `zb1105`
lives in a sibling repository: `../vector_flow_maps/`.

## Contents

```
miscellaneous_analyses/
├── README.md                 # this file
├── pyproject.toml            # Python deps (uv / pip)
├── requirements.txt          # pip fallback
├── environment.yml           # conda alternative
├── uv.lock                   # pinned Python lockfile
├── renv.lock                 # pinned R lockfile (renv)
├── .gitignore
│
├── methods/                  # manuscript methods write-up
│   ├── gastrulation_dynamics_methods.tex
│   └── gastrulation_dynamics_methods.pdf
│
├── scripts/                  # all analysis code, grouped by purpose
│   ├── viewers/              # napari-based interactive viewers (orientation, sphere fit, ROI)
│   │   ├── embryo_viewer.py          # TrackMate CSV input → oriented tracks
│   │   ├── ultrack_viewer.py         # ultrack-format CSV input
│   │   └── nuclei_tracks_viewer.py   # post-orientation viewer
│   ├── preprocessing/        # one-shot preprocessing
│   │   ├── nuclear_stats.py          # per-timepoint / per-slice nuclear stats from segmentation
│   │   └── export_to_mastodon.py     # Mastodon-format export
│   ├── dynamics/             # main gastrulation-dynamics analyses (per-species + cross-species)
│   │   ├── gastrulation_dynamics_medaka.R       # 10 figures per medaka recording
│   │   ├── gastrulation_dynamics_zebrafish.R    # 10 figures per zebrafish recording
│   │   └── gastrulation_dynamics_comparison.R   # 7 cross-species comparison figures
│   ├── oneshot_comparison/   # single-figure comparative panel (5 figures) + bulges/margins
│   │   ├── oneshot_comparison.R
│   │   ├── explore_medaka_bulge.R
│   │   └── explore_zeb_margin.R
│   ├── density/              # density / thickness / voronoi-cell-area analyses
│   │   ├── density_3timepoints.R                # 3-timepoint density (T-1h, T0, T+1h)
│   │   ├── medaka_density_analysis.R            # medaka density & flow figures
│   │   ├── zebrafish_density_analysis.R         # zebrafish density & flow figures
│   │   ├── fig_head_meso_vs_animalpole.R        # head-mesoderm vs animal-pole density
│   │   ├── density_marginbased_3tp.R            # voronoi-based margin density (came from the voronoi run)
│   │   └── voronoi_compute.py                   # voronoi tessellation pipeline (Python)
│   ├── surface_flow_maps/    # older Python surface-flow-map scripts (per-recording, NOT the Sep 2026 work)
│   │   ├── mk2508_surface_flow_maps.py
│   │   ├── mk2508_surface_flow_maps_1h_back.py
│   │   ├── mk2508_margin_density_1h_back.py
│   │   ├── mk3005_surface_flow_maps.py
│   │   ├── mk2805_surface_flow_maps_1h_back.py
│   │   └── mk2805_internalisation_crosssection.py
│   └── ultrack/
│       └── ultrack.def                  # ultrack configuration (detection + tracking)
│
├── data/
│   ├── raw/                  # raw input CSVs (TrackMate / ultrack exports)
│   │   ├── medaka_mk2508_spots.csv          # medaka 25-08-2025 (low-res) — TrackMate spots export
│   │   ├── medaka_mk2508_tracks.csv         # medaka 25-08-2025 (low-res) — TrackMate tracks export
│   │   └── medaka_mk2805_minimal_tracks.csv  # medaka 28-05-2025 minimal recording — tracks only
│   ├── medaka_mk3005_highres/        # medaka 30-05-2025 high-res raw tracks + ultrack artefacts
│   ├── medaka_mk2508_lowres/         # medaka 25-08-2025 low-res raw tracks + ultrack artefacts
│   ├── oriented_medaka_ultrack/      # oriented medaka tracks (output of viewer) — drives dynamics analyses
│   └── oriented_zebrafish_ultrack/   # oriented zebrafish tracks (output of viewer) — drives dynamics analyses
│
├── nuclear_stats/            # per-timepoint / per-slice / global nuclear stats (from nuclear_stats.py)
│   ├── medaka/
│   │   ├── nuclear_stats_global.tsv
│   │   ├── nuclear_stats_per_slice.tsv
│   │   ├── nuclear_stats_per_time.tsv
│   │   └── plots/
│   └── zebrafish/
│       ├── nuclear_stats_global.tsv
│       ├── nuclear_stats_per_slice.tsv
│       └── nuclear_stats_per_time.tsv
│
└── results/                  # all generated outputs (PDFs, CSVs, intermediate data)
    ├── medaka_dynamics/             # = analysis_output_medaka/  (10 dynamics PDFs + mk3005 surface flow maps)
    ├── zebrafish_dynamics/          # = analysis_output_zebrafish/  (10 dynamics PDFs)
    ├── zebrafish_0511/              # = analysis_output_zebrafish_05112025/  (older zebrafish dataset dynamics)
    ├── cross_species_comparison/    # = analysis_output_comparison/  (7 cross-species comparison PDFs)
    ├── oneshot_comparison/          # = analysis_output_oneshot/  (oneshot figures: margin density, flow, ingression vectors, matched tracks)
    ├── density_3timepoints/         # = analysis_output_density_3tp_20260625/  (3-tp density, Jun 25 2026)
    ├── density_marginbased/         # = analysis_output_density_marginbased_20260630/  (voronoi-based, Jun 30 2026)
    ├── medaka_mk2508/               # = analysis_output_medaka_25082025/  (mk2508 input data + surface flow maps + margin density)
    ├── medaka_mk2805_minimal/       # = analysis_output_2805_minimal/  (mk2805 minimal: surface flow maps + internalisation cross-section)
    └── zebrafish_oriented/          # = analysis_output/  (oriented zebrafish data, viewer output)
```

## Quick start

### Python environment

```bash
# Recommended (uv)
uv sync

# OR pip
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# OR conda
conda env create -f environment.yml
conda activate spim_medaka_zebrafish_trajectory
```

### R environment

```r
renv::restore()
```

## Pipeline

| Step | Script | Language | Purpose |
|------|--------|----------|---------|
| 1 | `scripts/viewers/embryo_viewer.py` | Python (napari) | Interactive 4D viewer: orientation, sphere fit, margin selection, export |
| 1′ | `scripts/viewers/ultrack_viewer.py` | Python (napari) | Variant of viewer for ultrack-format data |
| 2 | `scripts/preprocessing/nuclear_stats.py` | Python | Per-timepoint nuclear density, size, internuclear distance from segmentation labels |
| 3 | `scripts/dynamics/gastrulation_dynamics_medaka.R` | R | Medaka-specific dynamics analysis (10 figures) |
| 3′ | `scripts/dynamics/gastrulation_dynamics_zebrafish.R` | R | Zebrafish-specific dynamics analysis (10 figures) |
| 4 | `scripts/dynamics/gastrulation_dynamics_comparison.R` | R | Cross-species comparison (7 figures) |
| S1 | `scripts/density/medaka_density_analysis.R` | R | Supplementary: density & flow figures (medaka) |
| S2 | `scripts/density/zebrafish_density_analysis.R` | R | Supplementary: density & flow figures (zebrafish) |
| S3 | `scripts/density/density_3timepoints.R` | R | Supplementary: 3-timepoint density analysis |
| S4 | `scripts/density/density_marginbased_3tp.R` + `voronoi_compute.py` | R + Python | Supplementary: voronoi-based margin density |
| S5 | `scripts/oneshot_comparison/oneshot_comparison.R` | R | Supplementary: one-shot cross-species panel |
| O1 | `scripts/surface_flow_maps/*` | Python | Older per-recording surface flow maps (mk2508, mk3005, mk2805) |

## Running the napari 4D viewer (Step 1)

```bash
# TrackMate export (medaka mk2508)
uv run python scripts/viewers/embryo_viewer.py data/raw/medaka_mk2508_spots.csv \
    --tracks data/raw/medaka_mk2508_tracks.csv

# ultrack-format export
uv run python scripts/viewers/ultrack_viewer.py data/oriented_*_ultrack/oriented_tracks_*.csv
```

**Inside the viewer:**
1. Use the **time slider** to scrub through frames — tracks form progressively.
2. **Pick Animal Pole** → click a nucleus → **Pick Dorsal** → click a nucleus.
3. **Orient Embryo** → standardises axes (AP → +Y, Dorsal → +X).
4. **Fit Sphere** → cap-aware sphere fit for SPIM thin-slab data.
5. **Apply Margin** → flag spots within a latitude band (adjustable sliders).
6. **Colour buttons** → depth / latitude / frame / margin overlay.
7. **Export** → writes oriented CSVs + `sphere_params.csv` to `data/oriented_*_ultrack/`.

## Running the R analyses

```r
source("scripts/dynamics/gastrulation_dynamics_medaka.R")     # → results/medaka_dynamics/
source("scripts/dynamics/gastrulation_dynamics_zebrafish.R")  # → results/zebrafish_dynamics/
source("scripts/dynamics/gastrulation_dynamics_comparison.R") # → results/cross_species_comparison/
```

## Reproducibility

- **Python**: `uv sync` reproduces the exact environment from `pyproject.toml` + `uv.lock`.
- **R**: `renv::restore()` reproduces the R environment from `renv.lock`.
- **Data**: Raw CSVs in `data/raw/` are too large for git (see `.gitignore`).
  Store them on a shared drive / Zenodo / figshare.

## Provenance

This folder contains scripts, raw inputs, oriented data, and per-recording
outputs that were previously intermixed at the repository root. The
`archive/` folder at the workspace root retains the older scripts
(trackmate_*.R, orientation_interactive*.R, etc.) for historical context.
