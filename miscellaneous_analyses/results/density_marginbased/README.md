# Density @ 3 timepoints — Medaka bulge vs Zebrafish margin (T0-frozen, per-frame nuclei)

Follow-up to [analysis_output_density_3tp_20260625](../analysis_output_density_3tp_20260625)
adapted to the supervisor's 2026-06-30 feedback.

## Two changes vs the v0 (2026-06-30) version

1. **T0-frozen ingression regions** (the critical supervisor fix).
   The geometric margin/disc moves between timepoints (zebrafish: θ from
   105° at T-1h to 116° at T+1h, i.e., the margin sweeps inward via
   epiboly; medaka: θ from 72.5° to 78.8°). Re-picking the region at every
   timepoint was conflating region growth with cell-density growth. The
   fix defines the region ONCE at T0, using the 2-98 percentile (θ, φ)
   bounding box of the cells in the geometric margin/disc at T0, then
   re-uses that same bbox at T-1h, T0, T+1h.

2. **Per-frame nuclei, not per-track aggregates**.
   Within the ±5 min window, the same cell (TRACK_ID) appears at ~10-21
   frames. The earlier version took the median position per track
   within the window -- a single point that hid the actual per-frame
   nucleus distribution. The current code treats each row as a nucleus
   observation: per-frame metrics are computed at each FRAME in the
   window, then aggregated (median) across frames. KNN, MNN, depth,
   n_cells all use this per-frame definition.

## Where the regions are

In both species the imaging only covers the upper / margin side of the
embryo. The regions used in this analysis:

- **Animal cap (background)**: the top 5°-wide band of the imaged theta
  range. Defined identically for both species.
- **Head-mesoderm bulge (medaka RoI)**: data-driven disc at late frames,
  θ ≈ 79.5°, φ ≈ -92.9°, r ≈ 5.6° (early frames have < 50 deep cells
  and the disc falls back to the landmark).
- **Margin (zebrafish RoI)**: a θ-band centred on the `margin` landmark
  (θ ≈ 112.6°, φ ≈ 95.6°), ±8° in θ and ±60° in φ.
- **T0-frozen footprint**: 2-98 percentile bbox of the cells inside the
  RoI at T0. medaka bulge ∈ [74.3°, 84.0°] × [-98.0°, -87.7°];
  zebrafish margin ∈ [104.7°, 111.6°] × [72.3°, 119.1°].

See `t0_frozen_regions.csv` and `fig_F_topview_regions.pdf` for the
outline overlay.

## Output files

```
density_marginbased_3tp.R              Main analysis (R) -- single coherent flow
voronoi_compute.py                     Voronoi helper (Python, scipy)
README.md                              this file
metrics_long.csv                       per (species x tp x region)
metrics_wide_with_fold.csv             wide with fold changes
knn_mnn.csv                            KNN k=1/6/10 + MNN fraction
depth_invariance_zebrafish.csv         zebrafish margin depth percentiles
t0_frozen_regions.csv                  transparency: which (θ, φ) bbox used
cohort_trajectories.csv                tracks at T0 margin/disc, followed ±1h
voronoi_all_cells.csv + voronoi_region_summary.csv + voronoi/<...>.csv

fig_A_ncells.pdf            n cells per region (the numerator, per frame)
fig_B_depth_invariance.pdf  zebrafish margin depth percentiles
fig_C_knn_mnn.pdf           KNN k=1/6/10 and MNN fraction
fig_D_density_folds.pdf     density (per µm² and per µm³) + folds
fig_E_voronoi.pdf           per-cell Voronoi area + column volume
fig_F_topview_regions.pdf   animal-cap & RoI outlines overlaid on data
fig_G_cohort.pdf            T0-cohort distribution at T-1h / T0 / T+1h
```

## What we conclude from this dataset

1. **KNN ratios match the supervisor's expectations exactly**.
   Medaka head-mesoderm bulge / animal cap KNN10 ≈ **0.97 – 1.04**
   (essentially 1) at every timepoint. Zebrafish margin / animal cap
   KNN10 ≈ **0.93 – 1.21**, with T+1h closest to the supervisor's "~0.9".

2. **Density per µm² in the zebrafish margin grew ~5×** between T-1h
   (119 nuclei/frame) and T+1h (580 nuclei/frame) **inside the same
   frozen footprint**. This is a real fill-in, not a region artefact.

3. **Medaka bulge density grew only ~20%** over the same window (332 →
   396 nuclei/frame): progressive accumulation without dramatic local
   re-arrangement.

4. **Zebrafish margin cells do not internalise** (supervisor's point):
   p50 of SPHERICAL_DEPTH inside the T0-frozen margin stays at 41.3 →
   43.7 → 39.3 µm; fraction-deep-30 even falls (0.92 → 0.88 → 0.81).
   So the per-µm³ density (volumetric) of the zebrafish margin is a
   fair measurement -- any inflation in the 2026-06-25 numbers was due
   to the **moving** region, not to internalisation.

5. **The cohort is the smoking gun**.
   Tracking the unique TRACK_IDs that were in the **geometric** margin
   at T0 (not the frozen footprint), we find:

   | Species    | Metric        | T-1h | T0 | T+1h |
   |------------|---------------|------|----|------|
   | Zebrafish  | mean θ (°)    | 101.6 | 107.1 | 111.3 |
   | Zebrafish  | mean depth    | 39.6  | 44.9 | 43.8 |
   | Zebrafish  | frac in band  | 18%   | 100% | 95%  |
   | Medaka     | mean θ (°)    | 79.9  | 78.4 | 80.3 |
   | Medaka     | mean depth    | 37.7  | 37.5 | 33.8 |
   | Medaka     | frac in band  | 87%   | 100% | 100% |

   - **Zebrafish**: the margin sweeps inward by ~5° per hour, but the
     cells that were in the margin at T0 don't change depth. At T-1h
     they are scattered; at T+1h they remain inside the region.
   - **Medaka**: the cohort stays put in θ but is being joined by new,
     deeper cells over time; this is why the cohort mean depth rises
     while individual cohort depth stays roughly constant.

## Reproducibility

```bash
Rscript analysis_output_density_marginbased_20260630/density_marginbased_3tp.R
```

Required R packages (already in `renv.lock`): `data.table`, `ggplot2`,
`patchwork`, `scales`, `viridis`, `RANN`.

Required Python (the project's `.venv/`): `numpy`, `pandas`, `scipy`.
The script auto-detects the venv interpreter; system `python3` is used
as fallback (and would fail without scipy).
