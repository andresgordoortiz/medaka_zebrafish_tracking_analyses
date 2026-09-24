# spim_medaka_zebrafish_analyses

R project: trajectory and gastrulation-dynamics analyses of **SPIM (light-sheet) recordings of medaka and zebrafish embryos**.

. Key outputs are displacement/speed/straightness summaries, ingression diagnostics, nuclear density, orientation fields, and per-cell movement summaries around gastrulation landmarks (margin, YSL drift, internalisation, contraction).

## Layout

- **`vector_flow_maps/`** — subproject: binned vector-flow fields from oriented+filtered tracks for medaka `mk2508` and zebrafish `zb1105` (Sept 2026). Self-contained with its own `renv.lock`. See its README.
- **`miscellaneous_analyses/`** — unrelated Python side-analyses (own `environment.yml` / `pyproject.toml`).
- **`data/`, `exported_from_trackmate_gui/`, `high_res_mannualtracking/`, `output_medaka_28052025/`, `all_tracks/`, `nuclear_stats/`** — raw/processed inputs and TrackMate exports per dataset.

- **`archive/`** — superseded/working R scripts (TrackMate import, filtering, validation, orientation, density, gastrulation dynamics, Imaris comparison, manual-tracking reconciliation) plus prior `analysis_output_old/` snapshots. Used as the working scratch; canonical runs live in subprojects.

- **`environment.yml`, `requirements.txt`, `renv/`** — conda + pip + renv environments.

## Reproducing

- R analyses: `renv` is pinned; activate with `Rscript -e 'source("renv/activate.R")'`. The `.Rproj` file marks the project root.
- TrackMate / Imaris CSVs in `archive/all_tracks/` and `exported_from_trackmate_gui/` are the canonical inputs to `archive/trackmate_analysis.R` and downstream scripts.
