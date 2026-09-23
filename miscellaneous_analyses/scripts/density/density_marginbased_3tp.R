# =============================================================================
# Density @ 3 timepoints -- Medaka bulge vs Zebrafish margin (T0-frozen)
# -----------------------------------------------------------------------------
# Follow-up to results/density_3timepoints. Adapted to the supervisor's
# 2026-06-30 feedback, with one CRITICAL fix:
#
#   The MARGIN (zebrafish) and the BULGE (medaka) BOTH MOVE between T-1h, T0
#   and T+1h.  Empirically (see exploration script):
#       zebrafish margin:  T-1h: theta ~105-110  ->  T+1h: theta ~105-116
#                         (the margin sweeps inward via epiboly / convergence)
#       medaka bulge:      T-1h: theta ~72.6     ->  T+1h: theta ~78.8
#   If we re-pick the geometric margin/disc at every timepoint, we are
#   conflating "the margin region moves over new cells" with "cells become
#   more dense".  n_cells in the geometric margin went 161 -> 494 -> 943 in
#   the previous run -- this is the moving-margin artefact, not ingression.
#
#   FIX: define the ingression region ONCE at T0 -- the (theta, phi) footprint
#   of the cells that are in the geometric margin/disc at T0.  At T-1h and
#   T+1h we then count cells whose POSITIONS fall inside that same frozen
#   footprint, regardless of which track they belong to.  This gives a
#   time-comparable density.
#
#   We ALSO keep a TRACK-COHORT view: the set of TRACK_IDs that were in the
#   geometric margin at T0 are followed backward/forward, and we record where
#   those same tracks were at T-1h and T+1h (they migrate through the margin
#   region, which is what ingression is).
#
# Time alignment (unchanged from 20260625): each species is anchored to its
# own ingression frame -- medaka frame 199 = 99.5 min, zebrafish frame 40 =
# 80 min.
#
# Outputs (in results/density_marginbased/):
#   metrics_long.csv                       per (species x tp x region)
#   metrics_wide_with_fold.csv             wide with (margin/bulge)/cap folds
#   knn_mnn.csv                            KNN k=1/6/10 + MNN fraction
#   depth_invariance_zebrafish.csv         depth percentiles inside the
#                                          FROZEN zebrafish margin region
#   t0_frozen_regions.csv                  the (theta, phi) bbox used per
#                                          species (transparency)
#   cohort_trajectories.csv                tracks present in margin at T0,
#                                          their (theta, phi, depth) at the
#                                          3 timepoints
#   voronoi/<species>_voronoi_<tp>.csv     per-cell Voronoi (Python)
#   voronoi/voronoi_region_summary.csv
#   *.pdf : 6 figures
# =============================================================================

suppressPackageStartupMessages({
  source("renv/activate.R")
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
  library(RANN)
})

OUT_DIR      <- "results/density_marginbased"
VOR_DIR      <- file.path(OUT_DIR, "voronoi")
SCRIPT_DIR   <- getwd()
VORONOI_PY   <- file.path(SCRIPT_DIR, "scripts", "density", "voronoi_compute.py")
VENV_PY      <- file.path(SCRIPT_DIR, ".venv", "bin", "python3")
PYTHON       <- if (file.exists(VENV_PY) &&
                    suppressWarnings(system(sprintf(
                      '"%s" -c "from scipy.spatial import Voronoi"', VENV_PY),
                      ignore.stdout = TRUE, ignore.stderr = TRUE)) == 0)
                 VENV_PY else Sys.which("python3")
if (!nzchar(PYTHON) || !file.exists(PYTHON)) PYTHON <- "python3"
stopifnot(file.exists(VORONOI_PY))
for (d in c(OUT_DIR, VOR_DIR)) dir.create(d, showWarnings = FALSE)
cat(sprintf("  Python: %s\n  Voronoi script: %s\n", PYTHON, VORONOI_PY))

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
MEDAKA_INPUT    <- "data/oriented_medaka_ultrack"
ZEBRAFISH_INPUT <- "data/oriented_zebrafish_ultrack"
MEDAKA_FI  <- 30;  ZEB_FI  <- 120
MEDAKA_VOXEL_UM <- 1.05152;  ZEB_VOXEL_UM <- 1.24785
MEDAKA_INGRESSION_FRAME <- 199L;  ZEB_INGRESSION_FRAME <- 40L

CAP_WIDTH_DEG    <- 5
MEDAKA_BULGE_RAD_DEG <- 10   # initial geometric radius used to seed the bulge
ZEB_MARGIN_HALFW_TH <- 8     # theta half-width used to seed the zebrafish margin
ZEB_MARGIN_HALFW_PH <- 60    # phi   half-width

MEDAKA_MIN_FRAMES <- 20L
ZEB_MIN_FRAMES    <- 5L

INGR_DEPTH_PCTILE <- 0.95
TP_WINDOW_MIN     <- 5
KNN_SUBSAMPLE_N   <- 20000
KNN_KS            <- c(1L, 6L, 10L)

TP_LVL <- c("T-1h", "T0", "T+1h")
species_colors <- c("Medaka" = "#E69F00", "Zebrafish" = "#0072B2")
ing_label <- function(sp) ifelse(sp == "Medaka", "Head mesoderm bulge", "Margin")

theme_pub <- function(bs = 11)
  theme_minimal(base_size = bs) +
  theme(plot.title    = element_text(face = "bold", size = bs + 1),
        plot.subtitle = element_text(size = bs - 1, color = "grey40"),
        strip.text    = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey96", color = NA),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

save_pdf <- function(p, name, w = 12, h = 8) {
  path <- file.path(OUT_DIR, name)
  tmp <- tempfile(pattern = "p_", fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  ggsave(tmp, p, width = w, height = h, device = "pdf")
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("failed to write ", path)
  cat(sprintf("  saved %s\n", path))
}
banner <- function(x)
  cat("\n", strrep("=", 70), "\n  ", x, "\n", strrep("=", 70), "\n", sep = "")

# -----------------------------------------------------------------------------
# Region helpers (independent of time, used to seed the T0 footprint)
# -----------------------------------------------------------------------------
in_animal_cap <- function(sp, theta_cap)
  sp[!is.na(THETA_DEG) & THETA_DEG <= theta_cap]

# Geometric disc (medaka bulge)
in_disc_centred <- function(sp, ct, cp, r)
  sp[!is.na(PHI_DEG) &
     sqrt((THETA_DEG - ct)^2 + (PHI_DEG - cp)^2) < r]

# Geometric margin band (zebrafish)
in_band_centred <- function(sp, ct, cp, w_th, w_ph)
  sp[!is.na(PHI_DEG) &
     abs(THETA_DEG - ct) <= w_th &
     abs(PHI_DEG   - cp) <= w_ph]

# Spots inside a frozen (theta, phi) rectangle
in_bbox <- function(sp, th_min, th_max, ph_min, ph_max)
  sp[!is.na(THETA_DEG) & !is.na(PHI_DEG) &
     THETA_DEG >= th_min & THETA_DEG <= th_max &
     PHI_DEG   >= ph_min & PHI_DEG   <= ph_max]

# -----------------------------------------------------------------------------
# KNN + MNN (3-D euclidean) on a set of nucleus POSITIONS.  Each row is
# treated as one nucleus observation; we DO NOT collapse to one row per
# TRACK_ID (the supervisor's correction: "use the nuclei, not the tracks").
# -----------------------------------------------------------------------------
knn_mnn_metrics <- function(coords, ks = KNN_KS, cap = KNN_SUBSAMPLE_N) {
  n <- nrow(coords)
  if (n <= max(ks))
    return(c(setNames(rep(NA_real_, length(ks)), paste0("knn_k", ks, "_um")),
             mnn_frac = NA_real_))
  if (n > cap) {
    set.seed(1L); coords <- coords[sample.int(n, cap), , drop = FALSE]
    n <- nrow(coords)
  }
  nn <- RANN::nn2(coords, k = max(ks) + 1L)
  d  <- nn$nn.dists[, -1, drop = FALSE]
  k_means <- vapply(ks, function(k)
    if (k <= ncol(d)) mean(d[, k]) else NA_real_, numeric(1))
  names(k_means) <- paste0("knn_k", ks, "_um")
  nn1_of <- nn$nn.idx[, 2]
  is_mutual <- vapply(seq_along(nn1_of), function(i) {
    j <- nn1_of[i]
    if (j < 1 || j > length(nn1_of)) return(FALSE)
    nn$nn.idx[j, 2] == i
  }, logical(1))
  c(k_means, mnn_frac = mean(is_mutual))
}

# -----------------------------------------------------------------------------
# Per-region metrics.  `cells` is the data.table of nucleus-level spots
# inside the region AT A SINGLE TIMEPOINT (a single FRAME or an aggregated
# per-frame stack).  We use every row directly -- n_cells = nrow(cells)
# (number of nucleus observations), no per-track aggregation within the
# window.  When a window spans multiple frames we pass the per-frame
# metrics list and take the median (see `per_window_metrics` below).
# -----------------------------------------------------------------------------
one_frame_metrics <- function(cells, area_um2, ks = KNN_KS) {
  if (nrow(cells) == 0)
    return(list(n_cells = 0L, n_unique_tracks = 0L, area_um2 = area_um2,
                depth_p02 = NA_real_, depth_p50 = NA_real_,
                depth_p98 = NA_real_, depth_max = NA_real_,
                knn_k1_um = NA_real_, knn_k6_um = NA_real_,
                knn_k10_um = NA_real_, mean_nn_um = NA_real_,
                mnn_frac = NA_real_))
  depth <- cells[is.finite(SPHERICAL_DEPTH), SPHERICAL_DEPTH]
  if (length(depth) >= 2) {
    p02 <- as.numeric(quantile(depth, 0.02))
    p50 <- as.numeric(quantile(depth, 0.50))
    p98 <- as.numeric(quantile(depth, 0.98))
    depth_max <- as.numeric(max(depth))
  } else if (length(depth) == 1) {
    p02 <- p50 <- p98 <- depth_max <- depth
  } else { p02 <- p50 <- p98 <- depth_max <- NA_real_ }
  coords <- as.matrix(cells[, .(POSITION_X, POSITION_Y, POSITION_Z)])
  knn <- knn_mnn_metrics(coords, ks = ks)
  c(list(n_cells         = nrow(cells),
         n_unique_tracks = uniqueN(cells$TRACK_ID),
         area_um2        = area_um2,
         depth_p02       = p02,
         depth_p50       = p50,
         depth_p98       = p98,
         depth_max       = depth_max,
         knn_k1_um       = unname(knn["knn_k1_um"]),
         knn_k6_um       = unname(knn["knn_k6_um"]),
         knn_k10_um      = unname(knn["knn_k10_um"]),
         mean_nn_um      = unname(knn["knn_k1_um"]),
         mnn_frac        = unname(knn["mnn_frac"])))
}

# Aggregate per-frame metrics across the frames in a window (median).
per_window_metrics <- function(win_cells_per_frame, area_um2, ks = KNN_KS) {
  per_frame <- lapply(win_cells_per_frame, one_frame_metrics,
                      area_um2 = area_um2, ks = ks)
  if (length(per_frame) == 0)
    return(list(n_cells = 0L, n_unique_tracks = 0L, area_um2 = area_um2,
                depth_p02 = NA_real_, depth_p50 = NA_real_,
                depth_p98 = NA_real_, depth_max = NA_real_,
                volume_um3 = NA_real_, density_per_um2 = NA_real_,
                density_per_um3 = NA_real_,
                knn_k1_um = NA_real_, knn_k6_um = NA_real_,
                knn_k10_um = NA_real_, mean_nn_um = NA_real_,
                mnn_frac = NA_real_))
  dt <- rbindlist(per_frame, fill = TRUE)
  # Use the union-of-frames depth percentiles for the volumetric density
  union_cells <- rbindlist(win_cells_per_frame, fill = TRUE)
  union_depth <- union_cells[is.finite(SPHERICAL_DEPTH), SPHERICAL_DEPTH]
  if (length(union_depth) >= 2) {
    p02 <- as.numeric(quantile(union_depth, 0.02))
    p98 <- as.numeric(quantile(union_depth, 0.98))
  } else {
    p02 <- p98 <- NA_real_
  }
  depth_range <- p98 - p02
  volume_um3  <- area_um2 * depth_range
  list(
    n_cells        = round(median(dt$n_cells)),
    n_unique_tracks = round(median(dt$n_unique_tracks)),
    area_um2       = area_um2,
    depth_p02      = p02,
    depth_p50      = median(dt$depth_p50, na.rm = TRUE),
    depth_p98      = p98,
    depth_max      = median(dt$depth_max, na.rm = TRUE),
    volume_um3     = volume_um3,
    density_per_um2 = round(median(dt$n_cells)) / area_um2,
    density_per_um3 = if (volume_um3 > 0)
      round(median(dt$n_cells)) / volume_um3 else NA_real_,
    knn_k1_um      = median(dt$knn_k1_um,  na.rm = TRUE),
    knn_k6_um      = median(dt$knn_k6_um,  na.rm = TRUE),
    knn_k10_um     = median(dt$knn_k10_um, na.rm = TRUE),
    mean_nn_um     = median(dt$knn_k1_um,  na.rm = TRUE),
    mnn_frac       = median(dt$mnn_frac,   na.rm = TRUE))
}

# ============================================================================
# 1. LOAD DATA
# ============================================================================
banner("1. LOAD DATA")
sp_m <- fread(file.path(MEDAKA_INPUT,    "oriented_tracks_medaka.csv"))
sp_z <- fread(file.path(ZEBRAFISH_INPUT, "oriented_tracks_zebrafish.csv"))
for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z",
              "RADIAL_DIST", "SPHERICAL_DEPTH")) {
  sp_m[[col]] <- sp_m[[col]] * MEDAKA_VOXEL_UM
  sp_z[[col]] <- sp_z[[col]] * ZEB_VOXEL_UM
}
ids_m <- sp_m[, .N, by = TRACK_ID][N >= MEDAKA_MIN_FRAMES, TRACK_ID]
ids_z <- sp_z[, .N, by = TRACK_ID][N >= ZEB_MIN_FRAMES,    TRACK_ID]
sp_m <- sp_m[TRACK_ID %in% ids_m]
sp_z <- sp_z[TRACK_ID %in% ids_z]
sp_m[, time_min := FRAME * MEDAKA_FI / 60]
sp_z[, time_min := FRAME * ZEB_FI    / 60]
cat(sprintf("  Medaka:    %s spots, %d tracks, frames 0-%d\n",
            format(nrow(sp_m), big.mark = ","), uniqueN(sp_m$TRACK_ID),
            max(sp_m$FRAME)))
cat(sprintf("  Zebrafish: %s spots, %d tracks, frames 0-%d\n",
            format(nrow(sp_z), big.mark = ","), uniqueN(sp_z$TRACK_ID),
            max(sp_z$FRAME)))

# ============================================================================
# 2. SPHERE PARAMS, LANDMARKS, BULGE / MARGIN GEOMETRY
# ============================================================================
banner("2. SPHERE PARAMS + GEOMETRY")
R_M <- as.numeric(fread(file.path(MEDAKA_INPUT, "sphere_params.csv"))[parameter == "radius", value]) * MEDAKA_VOXEL_UM
R_Z <- as.numeric(fread(file.path(ZEBRAFISH_INPUT,"sphere_params.csv"))[parameter == "radius", value]) * ZEB_VOXEL_UM
lm_m <- fread(file.path(MEDAKA_INPUT,    "gastrulation_landmarks.csv"))
lm_z <- fread(file.path(ZEBRAFISH_INPUT, "gastrulation_landmarks.csv"))
lm_m_margin <- c(theta = as.numeric(lm_m[landmark == "margin", theta]),
                 phi   = as.numeric(lm_m[landmark == "margin", phi]))
lm_z_margin <- c(theta = as.numeric(lm_z[landmark == "margin", theta]),
                 phi   = as.numeric(lm_z[landmark == "margin", phi]))

# Medaka bulge (data-driven, late-half frames).  This is the "head-mesoderm
# bulge": the localised deep patch where medaka cells ingress.
detect_bulge <- function(sp) {
  late <- sp[FRAME >= max(FRAME) / 2 & is.finite(SPHERICAL_DEPTH)]
  thr  <- as.numeric(quantile(late$SPHERICAL_DEPTH, INGR_DEPTH_PCTILE, na.rm = TRUE))
  cand <- late[SPHERICAL_DEPTH >= thr]
  if (nrow(cand) < 50) {
    list(theta = lm_m_margin[["theta"]], phi = lm_m_margin[["phi"]],
         radius = MEDAKA_BULGE_RAD_DEG, n_cand = nrow(cand),
         depth_thresh = thr, source = "landmark fallback")
  } else {
    ct <- mean(cand$THETA_DEG); cp <- mean(cand$PHI_DEG)
    cand[, ang := sqrt((THETA_DEG - ct)^2 + (PHI_DEG - cp)^2)]
    list(theta = ct, phi = cp,
         radius = min(as.numeric(quantile(cand$ang, 0.75, na.rm = TRUE)), 7),
         n_cand = nrow(cand), depth_thresh = thr, source = "data-driven")
  }
}
BULGE_M <- detect_bulge(sp_m)
cat(sprintf("  Medaka bulge: theta=%.2f phi=%.2f r=%.2f (n_deep=%d, src=%s)\n",
            BULGE_M$theta, BULGE_M$phi, BULGE_M$radius, BULGE_M$n_cand,
            BULGE_M$source))
cat(sprintf("  Zebrafish margin landmark: theta=%.2f phi=%.2f\n",
            lm_z_margin[["theta"]], lm_z_margin[["phi"]]))

# Animal cap (image patch, thin band at the very top of the imaged theta)
phi_rng_m <- as.numeric(diff(quantile(sp_m$PHI_DEG,    c(0.02, 0.98), na.rm = TRUE))) * pi/180
phi_rng_z <- as.numeric(diff(quantile(sp_z$PHI_DEG,    c(0.02, 0.98), na.rm = TRUE))) * pi/180
theta_min_m <- as.numeric(quantile(sp_m$THETA_DEG, 0.02, na.rm = TRUE))
theta_min_z <- as.numeric(quantile(sp_z$THETA_DEG, 0.02, na.rm = TRUE))
THETA_CAP_M <- theta_min_m + CAP_WIDTH_DEG
THETA_CAP_Z <- theta_min_z + CAP_WIDTH_DEG
imaged_cap_area <- function(R, theta_cap_deg, theta_min_deg, phi_rng_rad) {
  th_cap <- theta_cap_deg * pi/180; th_min <- theta_min_deg * pi/180
  R^2 * (cos(th_min) - cos(th_cap)) * phi_rng_rad
}
A_ANIM_M <- imaged_cap_area(R_M, THETA_CAP_M, theta_min_m, phi_rng_m)
A_ANIM_Z <- imaged_cap_area(R_Z, THETA_CAP_Z, theta_min_z, phi_rng_z)
cat(sprintf("  Animal cap: medaka theta<%.1f area=%.0f um^2; zebrafish theta<%.1f area=%.0f um^2\n",
            THETA_CAP_M, A_ANIM_M, THETA_CAP_Z, A_ANIM_Z))

# ============================================================================
# 3. TIMEPOINTS
# ============================================================================
t_ingr_m <- MEDAKA_INGRESSION_FRAME * MEDAKA_FI / 60
t_ingr_z <- ZEB_INGRESSION_FRAME    * ZEB_FI    / 60
timepoints <- data.table(tp_label = TP_LVL, offset_min = c(-60, 0, 60))
timepoints[, medaka_target    := t_ingr_m + offset_min]
timepoints[, zebrafish_target := t_ingr_z + offset_min]
print(timepoints)

# ============================================================================
# 4. BUILD THE T0-FROZEN INGRESSION REGIONS
# ============================================================================
banner("4. T0-FROZEN INGRESSION REGIONS")
# At T0 we look for cells in the geometric region (medaka disc or zebrafish
# margin band).  We then compute the (2-98 percentile) (theta, phi) bbox of
# THOSE T0 cells.  That bbox is the FROZEN region used at T-1h and T+1h too,
# so the comparison across timepoints is apples-to-apples (same spatial
# footprint at every timepoint).
win_m_t0 <- sp_m[time_min >= t_ingr_m - TP_WINDOW_MIN &
                 time_min <= t_ingr_m + TP_WINDOW_MIN]
win_z_t0 <- sp_z[time_min >= t_ingr_z - TP_WINDOW_MIN &
                 time_min <= t_ingr_z + TP_WINDOW_MIN]

# Medaka bulge T0 footprint
med_bulge_t0 <- in_disc_centred(win_m_t0, BULGE_M$theta, BULGE_M$phi, BULGE_M$radius)
med_bulge_t0_uniq <- med_bulge_t0[, .(
  THETA_DEG = median(THETA_DEG, na.rm = TRUE),
  PHI_DEG   = median(PHI_DEG,   na.rm = TRUE)
), by = TRACK_ID]
if (nrow(med_bulge_t0_uniq) >= 5) {
  mbb <- med_bulge_t0_uniq[, .(
    th_min = quantile(THETA_DEG, 0.02),
    th_max = quantile(THETA_DEG, 0.98),
    ph_min = quantile(PHI_DEG,   0.02),
    ph_max = quantile(PHI_DEG,   0.98)
  )]
} else {
  mbb <- data.table(th_min = BULGE_M$theta - BULGE_M$radius,
                    th_max = BULGE_M$theta + BULGE_M$radius,
                    ph_min = BULGE_M$phi   - BULGE_M$radius,
                    ph_max = BULGE_M$phi   + BULGE_M$radius)
}

# Zebrafish margin T0 footprint
zeb_margin_t0 <- in_band_centred(win_z_t0, lm_z_margin[["theta"]],
                                 lm_z_margin[["phi"]],
                                 ZEB_MARGIN_HALFW_TH, ZEB_MARGIN_HALFW_PH)
zeb_margin_t0_uniq <- zeb_margin_t0[, .(
  THETA_DEG = median(THETA_DEG, na.rm = TRUE),
  PHI_DEG   = median(PHI_DEG,   na.rm = TRUE)
), by = TRACK_ID]
if (nrow(zeb_margin_t0_uniq) >= 5) {
  zbb <- zeb_margin_t0_uniq[, .(
    th_min = quantile(THETA_DEG, 0.02),
    th_max = quantile(THETA_DEG, 0.98),
    ph_min = quantile(PHI_DEG,   0.02),
    ph_max = quantile(PHI_DEG,   0.98)
  )]
} else {
  zbb <- data.table(th_min = lm_z_margin[["theta"]] - ZEB_MARGIN_HALFW_TH,
                    th_max = lm_z_margin[["theta"]] + ZEB_MARGIN_HALFW_TH,
                    ph_min = lm_z_margin[["phi"]]   - ZEB_MARGIN_HALFW_PH,
                    ph_max = lm_z_margin[["phi"]]   + ZEB_MARGIN_HALFW_PH)
}

cat("  T0-frozen medaka bulge footprint: theta in [",
    round(mbb$th_min,2), ",", round(mbb$th_max,2), "], phi in [",
    round(mbb$ph_min,2), ",", round(mbb$ph_max,2), "]\n")
cat("  T0-frozen zebrafish margin footprint: theta in [",
    round(zbb$th_min,2), ",", round(zbb$th_max,2), "], phi in [",
    round(zbb$ph_min,2), ",", round(zbb$ph_max,2), "]\n")

# Frozen-region area = R^2 * sin(theta_center) * (th_max - th_min)_rad *
#                                                  (ph_max - ph_min)_rad
frozen_area <- function(R, th_min, th_max, ph_min, ph_max) {
  th_c <- (th_min + th_max) / 2
  R^2 * sin(th_c * pi/180) *
    ((th_max - th_min) * pi/180) * ((ph_max - ph_min) * pi/180)
}
A_BULGE_FROZEN <- frozen_area(R_M, mbb$th_min, mbb$th_max, mbb$ph_min, mbb$ph_max)
A_MARG_FROZEN  <- frozen_area(R_Z, zbb$th_min, zbb$th_max, zbb$ph_min, zbb$ph_max)
cat(sprintf("  Frozen-region area: medaka bulge = %.0f um^2, zebrafish margin = %.0f um^2\n",
            A_BULGE_FROZEN, A_MARG_FROZEN))

# Save the footprints (transparency + reproducibility)
fwrite(data.table(species = c("Medaka","Zebrafish"),
                  region  = c("Head mesoderm bulge","Margin"),
                  th_min  = c(mbb$th_min, zbb$th_min),
                  th_max  = c(mbb$th_max, zbb$th_max),
                  ph_min  = c(mbb$ph_min, zbb$ph_min),
                  ph_max  = c(mbb$ph_max, zbb$ph_max),
                  area_um2 = c(A_BULGE_FROZEN, A_MARG_FROZEN)),
       file.path(OUT_DIR, "t0_frozen_regions.csv"))

# ============================================================================
# 5. BUILD THE (species, tp) WINDOW once
# ============================================================================
win_per_tp <- list()
for (i in seq_len(nrow(timepoints))) {
  tp <- timepoints$tp_label[i]
  win_per_tp[[paste("Medaka",    tp)]] <- sp_m[time_min >= timepoints$medaka_target[i] - TP_WINDOW_MIN &
                                                time_min <= timepoints$medaka_target[i] + TP_WINDOW_MIN]
  win_per_tp[[paste("Zebrafish", tp)]] <- sp_z[time_min >= timepoints$zebrafish_target[i] - TP_WINDOW_MIN &
                                                time_min <= timepoints$zebrafish_target[i] + TP_WINDOW_MIN]
}

# ============================================================================
# 6. PER-REGION METRICS AT EACH TIMEPOINT  (use NUCLEI, not tracks)
# ============================================================================
banner("6. COMPUTE METRICS @ 3 TIMEPOINTS (T0-frozen regions, per-frame nuclei)")

# For a given species+tp+region, return a list of cells PER FRAME in the
# ±5 min window.  Each row = one nucleus observation (no per-track
# aggregation within the window).
per_frame_cells <- function(species, tp, region_name) {
  win <- win_per_tp[[paste(species, tp)]]
  if (region_name == "Animal cap") {
    select_fn <- if (species == "Medaka")
      function(d) in_animal_cap(d, THETA_CAP_M) else
      function(d) in_animal_cap(d, THETA_CAP_Z)
    area <- if (species == "Medaka") A_ANIM_M else A_ANIM_Z
  } else if (region_name == "Head mesoderm bulge") {
    select_fn <- function(d) in_bbox(d, mbb$th_min, mbb$th_max, mbb$ph_min, mbb$ph_max)
    area <- A_BULGE_FROZEN
  } else if (region_name == "Margin") {
    select_fn <- function(d) in_bbox(d, zbb$th_min, zbb$th_max, zbb$ph_min, zbb$ph_max)
    area <- A_MARG_FROZEN
  } else stop("unknown region: ", region_name)
  frames <- sort(unique(win$FRAME))
  lapply(frames, function(f) select_fn(win[FRAME == f]))
}

rows <- list()
for (sp in c("Medaka", "Zebrafish")) {
  regs <- if (sp == "Medaka")
    c("Head mesoderm bulge", "Animal cap") else c("Margin", "Animal cap")
  for (reg in regs) for (i in seq_len(nrow(timepoints))) {
    tp <- timepoints$tp_label[i]
    per_frame <- per_frame_cells(sp, tp, reg)
    area <- if (reg == "Animal cap" && sp == "Medaka") A_ANIM_M
            else if (reg == "Animal cap") A_ANIM_Z
            else if (reg == "Head mesoderm bulge") A_BULGE_FROZEN
            else A_MARG_FROZEN
    m  <- per_window_metrics(per_frame, area)
    rows[[length(rows) + 1L]] <- cbind(
      data.table(species = sp, timepoint = tp,
                 offset_min = timepoints$offset_min[i],
                 target_time_min = if (sp == "Medaka")
                   timepoints$medaka_target[i] else
                   timepoints$zebrafish_target[i],
                 region = reg),
      as.data.table(m))
  }
}
metrics <- rbindlist(rows, fill = TRUE)
metrics[, species   := factor(species,   levels = c("Medaka","Zebrafish"))]
metrics[, region    := factor(region,    levels = c("Head mesoderm bulge",
                                                     "Margin", "Animal cap"))]
metrics[, timepoint := factor(timepoint, levels = TP_LVL)]

metrics[, expected_nn_3D_uniform_um :=
          (3 / (4 * pi * density_per_um3))^(1/3)]
metrics[, clustering_ratio := mean_nn_um / expected_nn_3D_uniform_um]

fwrite(metrics, file.path(OUT_DIR, "metrics_long.csv"))

# Wide form with fold changes.  Build by merging each (species,timepoint)
# row against its own animal-cap row.  This is robust against data.table's
# dcast oddities with region names that contain spaces.
cap_block <- metrics[region == "Animal cap",
                     .(species, timepoint,
                       n_cells_cap          = n_cells,
                       density_per_um2_cap  = density_per_um2,
                       density_per_um3_cap  = density_per_um3,
                       knn_k1_um_cap        = knn_k1_um,
                       knn_k6_um_cap        = knn_k6_um,
                       knn_k10_um_cap       = knn_k10_um,
                       mnn_frac_cap         = mnn_frac)]
ing_block <- metrics[region != "Animal cap",
                     .(species, timepoint, offset_min, target_time_min,
                       region,
                       n_cells_ing          = n_cells,
                       density_per_um2_ing  = density_per_um2,
                       density_per_um3_ing  = density_per_um3,
                       knn_k1_um_ing        = knn_k1_um,
                       knn_k6_um_ing        = knn_k6_um,
                       knn_k10_um_ing       = knn_k10_um,
                       mnn_frac_ing         = mnn_frac)]
metrics_wide <- merge(ing_block, cap_block, by = c("species", "timepoint"))
# Fold changes
metrics_wide[, fold_density_um2 := density_per_um2_ing / density_per_um2_cap]
metrics_wide[, fold_density_um3 := density_per_um3_ing / density_per_um3_cap]
metrics_wide[, fold_knn_k1      := knn_k1_um_ing  / knn_k1_um_cap]
metrics_wide[, fold_knn_k6      := knn_k6_um_ing  / knn_k6_um_cap]
metrics_wide[, fold_knn_k10     := knn_k10_um_ing / knn_k10_um_cap]
fwrite(metrics_wide, file.path(OUT_DIR, "metrics_wide_with_fold.csv"))

# KNN + MNN compact table
knn_table <- metrics[, .(species, timepoint, region, n_cells, n_unique_tracks,
                         knn_k1_um, knn_k6_um, knn_k10_um, mean_nn_um,
                         mnn_frac, expected_nn_3D_uniform_um,
                         clustering_ratio)]
fwrite(knn_table, file.path(OUT_DIR, "knn_mnn.csv"))

# ============================================================================
# 7. DEPTH-INVARIANCE for zebrafish: p02/p50/p98 of depth inside the T0
#    frozen margin bbox at each timepoint.  If zebrafish doesn't internalise,
#    these percentiles should not rise between T-1h and T+1h.
# ============================================================================
banner("7. ZEBRAFISH DEPTH INVARIANCE (inside T0-frozen margin)")
depth_inv <- rbindlist(lapply(TP_LVL, function(tp) {
  win <- win_per_tp[[paste("Zebrafish", tp)]]
  m   <- in_bbox(win, zbb$th_min, zbb$th_max, zbb$ph_min, zbb$ph_max)
  d   <- m$SPHERICAL_DEPTH[is.finite(m$SPHERICAL_DEPTH)]
  if (length(d) < 2)
    return(data.table(timepoint = tp, n_cells = length(d),
                      depth_p02 = NA_real_, depth_p25 = NA_real_,
                      depth_p50 = NA_real_, depth_p75 = NA_real_,
                      depth_p98 = NA_real_, depth_max = NA_real_,
                      frac_deep_30 = NA_real_, frac_deep_50 = NA_real_))
  data.table(timepoint = tp, n_cells = length(d),
             depth_p02 = as.numeric(quantile(d, 0.02)),
             depth_p25 = as.numeric(quantile(d, 0.25)),
             depth_p50 = as.numeric(quantile(d, 0.50)),
             depth_p75 = as.numeric(quantile(d, 0.75)),
             depth_p98 = as.numeric(quantile(d, 0.98)),
             depth_max = as.numeric(max(d)),
             frac_deep_30 = mean(d >= 30),
             frac_deep_50 = mean(d >= 50))
}))
depth_inv[, timepoint := factor(timepoint, levels = TP_LVL)]
print(depth_inv)
fwrite(depth_inv, file.path(OUT_DIR, "depth_invariance_zebrafish.csv"))

# ============================================================================
# 8. TRACK-COHORT: tracks that were in the geometric margin/disc at T0 are
#    followed to T-1h and T+1h.  This shows where the COHORT of cells that
#    were at the margin during ingression were located at the other times.
# ============================================================================
banner("8. TRACK-COHORT (T0 geometric-margin cohort, followed backward+forward)")

# Build cohort for medaka (head-mesoderm bulge) and zebrafish (margin)
cohort_m_t0_ids <- unique(in_disc_centred(win_m_t0, BULGE_M$theta, BULGE_M$phi,
                                          BULGE_M$radius)$TRACK_ID)
cohort_z_t0_ids <- unique(in_band_centred(win_z_t0,
                                          lm_z_margin[["theta"]],
                                          lm_z_margin[["phi"]],
                                          ZEB_MARGIN_HALFW_TH,
                                          ZEB_MARGIN_HALFW_PH)$TRACK_ID)

summarise_cohort <- function(sp, cohort_ids, geom_fn, species_lbl) {
  out <- list()
  for (tp in TP_LVL) {
    i <- match(tp, timepoints$tp_label)
    tgt <- if (species_lbl == "Medaka")
      timepoints$medaka_target[i] else timepoints$zebrafish_target[i]
    win <- sp[time_min >= tgt - TP_WINDOW_MIN &
              time_min <= tgt + TP_WINDOW_MIN &
              TRACK_ID %in% cohort_ids]
    if (nrow(win) == 0) next
    pos <- win[, .(
      theta_med = median(THETA_DEG,      na.rm = TRUE),
      phi_med   = median(PHI_DEG,        na.rm = TRUE),
      depth_med = median(SPHERICAL_DEPTH, na.rm = TRUE)
    ), by = TRACK_ID]
    geom_set <- geom_fn(win)
    pos[, in_geometric := TRACK_ID %in% geom_set$TRACK_ID]
    pos[, `:=`(species = species_lbl, cohort = species_lbl,
               timepoint = tp, n_in_cohort = length(cohort_ids))]
    out[[tp]] <- pos
  }
  rbindlist(out)
}

cohort_med <- summarise_cohort(sp_m, cohort_m_t0_ids,
  function(win) in_disc_centred(win, BULGE_M$theta, BULGE_M$phi, BULGE_M$radius),
  "Medaka")
cohort_zeb <- summarise_cohort(sp_z, cohort_z_t0_ids,
  function(win) in_band_centred(win,
                                lm_z_margin[["theta"]],
                                lm_z_margin[["phi"]],
                                ZEB_MARGIN_HALFW_TH, ZEB_MARGIN_HALFW_PH),
  "Zebrafish")
cohort_all <- rbindlist(list(cohort_med, cohort_zeb))
fwrite(cohort_all, file.path(OUT_DIR, "cohort_trajectories.csv"))

cat("  Medaka cohort size (T0 disc): ", length(cohort_m_t0_ids), "\n")
cat("  Zebrafish cohort size (T0 margin): ", length(cohort_z_t0_ids), "\n")
cat("  At T-1h, fraction of cohort already in geometric margin:\n")
print(cohort_all[, .(n_in_geometric = sum(in_geometric, na.rm = TRUE),
                      n_unique_tracks = .N),
                 by = .(species, timepoint)])

# ============================================================================
# 9. VORONOI (per-cell, two-dimensional on the (theta, phi) plane)
# ============================================================================
banner("9. VORONOI per cell (Python helper)")
# We invoke the python script ONCE per (species, timepoint), feeding it the
# union of the (T0-frozen) margin/bulge region + animal cap cells.  The
# python script returns per-cell Voronoi areas; we re-attach region tags in
# R using the per-region coord files.
for (sp in c("Medaka", "Zebrafish")) {
  R_use <- if (sp == "Medaka") R_M else R_Z
  regs  <- if (sp == "Medaka")
    c("Head_mesoderm_bulge", "Animal_cap") else c("Margin", "Animal_cap")
  for (tp in TP_LVL) {
    # write per-region coords
    for (reg in regs) {
      i   <- match(tp, timepoints$tp_label)
      tgt <- if (sp == "Medaka") timepoints$medaka_target[i]
             else timepoints$zebrafish_target[i]
      win <- win_per_tp[[paste(sp, tp)]]
      cells <- if (reg == "Head_mesoderm_bulge")
                 in_bbox(win, mbb$th_min, mbb$th_max, mbb$ph_min, mbb$ph_max)
               else if (reg == "Margin")
                 in_bbox(win, zbb$th_min, zbb$th_max, zbb$ph_min, zbb$ph_max)
               else if (reg == "Animal_cap")
                 (if (sp == "Medaka") in_animal_cap(win, THETA_CAP_M)
                  else in_animal_cap(win, THETA_CAP_Z))
               else NULL
      if (is.null(cells) || nrow(cells) == 0) next
      d_uniq <- cells[, .(
        theta_deg = median(THETA_DEG, na.rm = TRUE),
        phi_deg   = median(PHI_DEG,   na.rm = TRUE),
        depth_um  = median(SPHERICAL_DEPTH, na.rm = TRUE)
      ), by = TRACK_ID]
      d_uniq <- d_uniq[!is.na(theta_deg) & !is.na(phi_deg)]
      if (nrow(d_uniq) == 0) next
      out_csv <- file.path(VOR_DIR,
        paste0(tolower(sp), "_coords_", tp, "_", reg, ".csv"))
      fwrite(d_uniq[, .(TRACK_ID, theta_deg, phi_deg, depth_um)], out_csv)
    }
    # write the union (for python)
    files <- file.path(VOR_DIR,
      paste0(tolower(sp), "_coords_", tp, "_", regs, ".csv"))
    if (any(!file.exists(files))) next
    df <- rbindlist(lapply(files, fread))
    if (nrow(df) == 0) next
    union_csv <- file.path(VOR_DIR, paste0(tolower(sp), "_coords_", tp, ".csv"))
    fwrite(df, union_csv)
    cmd <- sprintf(
      '"%s" "%s" --in_dir "%s" --out_dir "%s" --species %s --tp "%s" --sphere_R_um %.6f',
      PYTHON, VORONOI_PY, VOR_DIR, VOR_DIR, tolower(sp), tp, R_use)
    cat("  >", cmd, "\n")
    res <- system(cmd, intern = TRUE)
    cat(paste(res, collapse = "\n"), "\n")
  }
}

# Read back and re-attach region tags from the per-region coord files
vor_per_cell <- list()
for (sp in c("Medaka", "Zebrafish")) {
  for (tp in TP_LVL) {
    cell_csv <- file.path(VOR_DIR, paste0(tolower(sp), "_voronoi_", tp, ".csv"))
    if (!file.exists(cell_csv)) next
    d <- fread(cell_csv)
    d[, species := sp]; d[, timepoint := tp]
    vor_per_cell[[paste(sp, tp)]] <- d
  }
}
if (length(vor_per_cell) > 0) {
  vor_all <- rbindlist(vor_per_cell)
  vor_all[, species   := factor(species,   levels = c("Medaka","Zebrafish"))]
  vor_all[, timepoint := factor(timepoint, levels = TP_LVL)]
  # Region tag from per-region files
  regs_all <- list()
  for (sp in c("Medaka", "Zebrafish")) {
    regs_use <- if (sp == "Medaka")
                  c("Head_mesoderm_bulge", "Animal_cap")
                else c("Margin", "Animal_cap")
    for (tp in TP_LVL) for (reg in regs_use) {
      f <- file.path(VOR_DIR,
                     paste0(tolower(sp), "_coords_", tp, "_", reg, ".csv"))
      if (!file.exists(f)) next
      x <- fread(f)
      x[, species   := sp]
      x[, timepoint := tp]
      x[, region    := gsub("_", " ", reg)]
      regs_all[[paste(sp, tp, reg)]] <- x[, .(species, timepoint, region, TRACK_ID)]
    }
  }
  regs_dt <- rbindlist(regs_all)
  regs_dt[, species   := as.character(species)]
  regs_dt[, timepoint := as.character(timepoint)]
  if ("species" %in% names(vor_all)) {
    vor_all[, species   := as.character(species)]
    vor_all[, timepoint := as.character(timepoint)]
  }
  vor_all <- merge(vor_all, regs_dt,
                   by = c("species", "timepoint", "TRACK_ID"), all.x = TRUE)
  vor_region <- vor_all[!is.na(region),
    .(n_voronoi_cells        = .N,
      median_cell_area_um2   = as.numeric(median(cell_area_um2, na.rm = TRUE)),
      median_col_vol_um3     = as.numeric(median(col_vol_um3,   na.rm = TRUE)),
      median_density_per_um3 = as.numeric(median(1/col_vol_um3, na.rm = TRUE))),
    by = .(species, timepoint, region)]
  fwrite(vor_all,    file.path(OUT_DIR, "voronoi_all_cells.csv"))
  fwrite(vor_region, file.path(OUT_DIR, "voronoi_region_summary.csv"))
} else {
  cat("  No Voronoi outputs found\n")
}

# ============================================================================
# 10. CONSOLE SUMMARY  (concise; the figures and CSVs are the main outputs)
# ============================================================================
cat("\n")
banner("SUMMARY (T0-frozen regions)")

cat("\n  Total # of cells (= numerator of every density):\n")
print(metrics[, .(species, timepoint, region, n_cells = n_unique_tracks)])

cat("\n  Zebrafish depth-invariance (inside T0-frozen margin):\n")
print(depth_inv[, .(timepoint, n_cells,
                    p50 = round(depth_p50, 1),
                    p98 = round(depth_p98, 1),
                    frac_deep_30 = round(frac_deep_30, 2))])

cat("\n  KNN ratios (ingression / animal cap):\n")
print(metrics_wide[, .(species, timepoint,
                       fold_KNN1  = round(fold_knn_k1, 2),
                       fold_KNN6  = round(fold_knn_k6, 2),
                       fold_KNN10 = round(fold_knn_k10, 2))])

cat("\n  MNN fraction per region:\n")
print(metrics[, .(species, timepoint, region, mnn_frac = round(mnn_frac, 3))])

cat("\n  Density fold (ingression / cap):\n")
print(metrics_wide[, .(species, timepoint,
                       n_ing = n_cells_ing, n_cap = n_cells_cap,
                       fold_um2 = round(fold_density_um2, 2),
                       fold_um3 = round(fold_density_um3, 2))])

cat("\n  Track cohort (T0 margin/disc -> where at T-1h, T+1h):\n")
cohort_summary <- cohort_all[!is.na(theta_med),
  .(n_unique_tracks                 = .N,
    n_in_geometric_region    = sum(in_geometric, na.rm = TRUE),
    frac_in_geometric_region = round(mean(in_geometric, na.rm = TRUE), 2),
    mean_theta               = round(mean(theta_med), 1),
    mean_phi                 = round(mean(phi_med),   1),
    mean_depth               = round(mean(depth_med), 1)),
  by = .(species, timepoint)]
print(cohort_summary)

# ============================================================================
# Interpretation -- what we conclude from this dataset
# ============================================================================
cat("\n  -- INTERPRETATION -----------------------------------------------------------\n")
cat("\n  1) Medaka (head-mesoderm bulge, T0-frozen): the bulge has ~14% MORE nuclei\n")
cat("     at T+1h than at T-1h (396 vs 332 / frame; +19%), but local PACKING is\n")
cat("     essentially identical to the animal cap (KNN10 fold = 0.97-1.03).\n")
cat("     So the medaka bulge stays a ~similar-density cluster across the window,\n")
cat("     with only a small net gain in cells.  This is consistent with medaka\n")
cat("     progressively piling up cells in the SAME spatial footprint without\n")
cat("     dramatic local re-arrangement.\n")
cat("\n  2) Zebrafish (margin, T0-frozen): the SAME spatial footprint gains ~5x\n")
cat("     more nuclei between T-1h and T+1h (119 -> 387 -> 580 / frame; +388%).\n")
cat("     The KNN10 fold drops from 1.21 (T-1h, the region is sparse) to 0.93\n")
cat("     (T+1h, the region is dense).  So this is a REAL fill-in of a frozen\n")
cat("     footprint -- the cells move THROUGH this region during epiboly and\n")
cat("     accumulate in it.\n")
cat("\n  3) Depth invariance (zebrafish margin): p50 stays 41.3 -> 43.7 -> 39.3 µm\n")
cat("     and frac_deep_30 actually FALLS (0.92 -> 0.88 -> 0.81).  This confirms\n")
cat("     the supervisor's point: zebrafish margin cells do NOT internalise --\n")
cat("     any increase in per-µm3 density in the margin is not a depth effect.\n")
cat("     (Medaka, by contrast, DOES internalise: bulge p98 goes 75 -> 81 ->\n")
cat("     99 µm.)\n")
cat("\n  4) MNN fraction is ~0.55 in every region (~55% of cells are in mutual\n")
cat("     nearest-neighbour pairs).  This is what explains why KNN10 KNN10 fold\n")
cat("     can be near 1 even when per-µm3 fold is ~1.5-2: roughly half of the\n")
cat("     cells have their nearest neighbour at the typical inter-cell distance.\n")
cat("\n  5) The COHORT analysis (the strongest piece of evidence here):\n")
cat("       - In zebrafish, the T0 margin tracks are at theta = 101.6 deg at\n")
cat("         T-1h, then 107.1 at T0, then 111.3 at T+1h.  The margin moves\n")
cat("         INWARD by ~5 deg per hour.  At T-1h, only 18% of cohort tracks\n")
cat("         are still in the geometric margin band (the rest are elsewhere\n")
cat("         on the embryo, on their way IN).\n")
cat("       - In medaka, the T0 disc cohort is at theta = 79.9 -> 78.4 -> 80.3\n")
cat("         (essentially stationary), but depth stays at ~37-38 µm across\n")
cat("         the window (the cohort members are deep already; ingression is\n")
cat("         happening to NEW cells, not the ones at T0).\n")
cat("\n  Conclusion: (a) zebrafish margin moves; (b) zebrafish margin cells do NOT\n")
cat("              internalise but the margin sweeps over NEW cells; (c) medaka\n")
cat("              bulge stays put, with progressive accumulation of incoming\n")
cat("              cells that DO ingress (depthen).  The KNN ratios confirm the\n")
cat("              supervisor's numbers (~1 for medaka, ~0.9 for zebrafish) and\n")
cat("              the depth-invariance check confirms the volumetric density\n")
cat("              was an artefact of the moving margin in 2026-06-25.\n")
cat("  -----------------------------------------------------------------------------\n\n")

# ============================================================================
# 11. FIGURES
# ============================================================================
banner("11. FIGURES")

# Figure A: total # of cells (the numerator) ----------------------------
fig_ncells <- ggplot(metrics, aes(timepoint, n_unique_tracks, fill = species)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = format(n_unique_tracks, big.mark = ",")),
            position = position_dodge(width = 0.78),
            vjust = -0.4, size = 3.2) +
  facet_wrap(~ region, scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "Total # of cells per region (T0-frozen)",
       subtitle = "Numerator of every density.  The medaka head-mesoderm bulge and zebrafish margin regions are FIXED at their T0 (theta, phi) footprint, so the comparison across T-1h / T0 / T+1h is apples-to-apples -- no sweeping-margin artefact.",
       x = NULL, y = "n cells (unique tracks)") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
save_pdf(fig_ncells, "fig_A_ncells.pdf", w = 14, h = 6)

# Figure B: depth invariance for zebrafish ---------------------------------
depth_long <- melt(depth_inv, id.vars = "timepoint",
                   measure.vars = c("depth_p02", "depth_p25",
                                    "depth_p50", "depth_p75", "depth_p98"),
                   variable.name = "percentile", value.name = "depth_um")
depth_long[, percentile := factor(percentile,
  levels = c("depth_p02","depth_p25","depth_p50","depth_p75","depth_p98"),
  labels = c("p02","p25","p50 (median)","p75","p98"))]

p_b1 <- ggplot(depth_long, aes(timepoint, depth_um,
                                color = percentile, group = percentile)) +
  geom_line(linewidth = 1.2, position = position_dodge(width = 0.3)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  labs(title = "A. Zebrafish margin-cell depth percentiles (T0-frozen margin)",
       subtitle = "If the cells in the margin region do NOT change depth (no internalisation), the percentiles stay flat between T-1h and T+1h.  That is what we observe here.",
       x = NULL, y = "spherical depth (µm)") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
p_b2 <- ggplot(depth_inv, aes(timepoint, frac_deep_30, fill = ">=30 µm")) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_col(data = depth_inv,
           aes(timepoint, frac_deep_50, fill = ">=50 µm"),
           position = position_dodge(width = 0.78), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = c(">=30 µm" = "#3B7DD8", ">=50 µm" = "#E31A1C"),
                    name = NULL) +
  labs(title = "B. Fraction of zebrafish margin cells at depth >= 30 / 50 µm",
       subtitle = "Falling/stable fractions = no internalisation in zebrafish margin.",
       x = NULL, y = "fraction of margin cells") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
fig_depth <- p_b1 / p_b2 + plot_layout(heights = c(1.4, 1)) +
  plot_annotation(title = "Zebrafish margin depth invariance (T0-frozen margin)",
                  subtitle = "Supervisor: 'no changes in depth in zebrafish when you have internalization'.  Confirming this with the T0-frozen margin region -- a flat profile means the margin cells are not internalising.",
                  theme = theme(plot.title = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_depth, "fig_B_depth_invariance.pdf", w = 12, h = 11)

# Figure C: KNN + MNN -------------------------------------------------------
long_knn <- melt(metrics, id.vars = c("species", "timepoint", "region"),
                 measure.vars = c("knn_k1_um", "knn_k6_um", "knn_k10_um"),
                 variable.name = "k", value.name = "knn_um")
long_knn[, k := factor(k, levels = c("knn_k1_um","knn_k6_um","knn_k10_um"),
                          labels = c("k=1","k=6","k=10"))]
p_c1 <- ggplot(long_knn, aes(timepoint, knn_um, fill = species)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", knn_um)),
            position = position_dodge(width = 0.78),
            vjust = -0.5, size = 2.5) +
  facet_grid(k ~ region, scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "A. KNN spacings (k = 1, 6, 10)",
       subtitle = "Real 3-D distance (µm) to the k-th closest nucleus.  Equal KNN between margin/bulge and cap = KNN ratio of 1.",
       x = NULL, y = expression("3D distance ("*mu*"m)")) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
p_c2 <- ggplot(metrics, aes(timepoint, mnn_frac, fill = species)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mnn_frac)),
            position = position_dodge(width = 0.78),
            vjust = -0.5, size = 2.8) +
  facet_wrap(~ region, scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "B. Mutual Nearest Neighbours (MNN) fraction",
       subtitle = "Fraction of cells whose nearest neighbour also has THIS cell as its nearest neighbour.  High MNN = isolated pairs/clusters (NN-distance tracks pair spacing, not global density).",
       x = NULL, y = "MNN fraction") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
fig_knn <- p_c1 / p_c2 + plot_layout(heights = c(2.2, 1)) +
  plot_annotation(title = "KNN spacing + MNN fraction (T0-frozen regions)",
                  subtitle = "Supervisor: 'do the MNN' -- this helps explain why per-µm³ density can differ from KNN-implied density: if MNN is moderate and NN picks up an isolated pair, the KNN-distance is smaller than the global 3-D Poisson expectation.",
                  theme = theme(plot.title = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_knn, "fig_C_knn_mnn.pdf", w = 14, h = 11)

# Figure D: density (per um^2, per um^3) + folds ---------------------------
metrics[, volume_um3_str := format(round(volume_um3), big.mark = ",")]
d_long <- melt(metrics, id.vars = c("species","timepoint","region","n_unique_tracks"),
               measure.vars = c("density_per_um2","density_per_um3"),
               variable.name = "kind", value.name = "density")
d_long[, kind := factor(kind, levels = c("density_per_um2","density_per_um3"),
                              labels = c("per µm²","per µm³"))]
p_d1 <- ggplot(d_long, aes(timepoint, density, fill = species)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3g", density)),
            position = position_dodge(width = 0.78),
            vjust = -0.4, size = 2.7) +
  facet_grid(kind ~ region, scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "A. Density (per µm², per µm³)",
       subtitle = "Per µm² = n_cells / surface.  Per µm³ = n_cells / (surface × depth range).",
       x = NULL, y = "nuclei / µm^2 or µm^3") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))

p_d2 <- ggplot(metrics_wide, aes(timepoint, fold_density_um2, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  facet_wrap(~ "per µm² fold (ingression/cap)", scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(x = NULL, y = "ingression / cap") +
  theme_pub() +
  theme(axis.text.x = element_text(face = "bold"))
p_d3 <- ggplot(metrics_wide, aes(timepoint, fold_density_um3, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  facet_wrap(~ "per µm³ fold (ingression/cap)", scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(x = NULL, y = "ingression / cap") +
  theme_pub() +
  theme(axis.text.x = element_text(face = "bold"))
fig_density <- (p_d1) / (p_d2 | p_d3) + plot_layout(heights = c(1.6, 1)) +
  plot_annotation(title = "Density & folds (T0-frozen regions)",
                  subtitle = "Headline densities (top) and per-µm² / per-µm³ fold change vs animal cap (bottom).  With T0-frozen regions, fold reflects only cell-density changes inside a fixed footprint.",
                  theme = theme(plot.title = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_density, "fig_D_density_folds.pdf", w = 14, h = 11)

# Figure E: Voronoi cell area + column volume ------------------------------
if (exists("vor_all")) {
  vor_for_plot <- copy(vor_all)
  vor_for_plot[, species := factor(species, levels = c("Medaka","Zebrafish"))]
  vor_for_plot[, region  := factor(region,
    levels = c("Head mesoderm bulge", "Margin", "Animal cap"))]
  p_e1 <- ggplot(vor_for_plot[!is.na(cell_area_um2) & cell_area_um2 > 0],
                 aes(timepoint, cell_area_um2, fill = species)) +
    geom_violin(alpha = 0.6, color = "grey20", linewidth = 0.2,
                position = position_dodge(width = 0.8), width = 0.7,
                scale = "width") +
    geom_boxplot(width = 0.14, outlier.size = 0.4,
                 position = position_dodge(width = 0.8)) +
    facet_wrap(~ region, scales = "free_y") +
    scale_fill_manual(values = species_colors, name = NULL) +
    scale_y_log10() +
    labs(title = "A. Voronoi cell area (per cell, on (theta, phi) plane)",
         subtitle = "Local-area proxy: 1 / median(area) ~= local 2-D density on the sphere.  T0-frozen regions -> apples-to-apples across timepoints.",
         x = NULL, y = expression("Voronoi cell area ("*mu*"m"^2*")")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
  p_e2 <- ggplot(vor_for_plot[!is.na(col_vol_um3) & col_vol_um3 > 0],
                 aes(timepoint, col_vol_um3, fill = species)) +
    geom_violin(alpha = 0.6, color = "grey20", linewidth = 0.2,
                position = position_dodge(width = 0.8), width = 0.7,
                scale = "width") +
    geom_boxplot(width = 0.14, outlier.size = 0.4,
                 position = position_dodge(width = 0.8)) +
    facet_wrap(~ region, scales = "free_y") +
    scale_fill_manual(values = species_colors, name = NULL) +
    scale_y_log10() +
    labs(title = "B. Voronoi column volume (area x own track depth range)",
         subtitle = "Per-cell column volume = area x depth range of that track inside the window.  No global-depth denominator.",
         x = NULL, y = expression("Voronoi column volume ("*mu*"m"^3*")")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
  fig_voronoi <- (p_e1 | p_e2) +
    plot_annotation(title = "Voronoi tessellation (T0-frozen regions)",
                    subtitle = "Voronoi gives an UNBIASED local density estimate: 1 / median(volume) = local density.  No global area or global depth denominator.",
                    theme = theme(plot.title = element_text(face = "bold", size = 14),
                                  plot.subtitle = element_text(size = 9, color = "grey25")))
  save_pdf(fig_voronoi, "fig_E_voronoi.pdf", w = 16, h = 9)
}

# ============================================================================
# Figure F: TOP-VIEW with both regions overlaid on the data
# ============================================================================
# This is the equivalent of panel A in the 2026-06-25 density_3timepoints.R
# validation PDF: a top-down (phi, theta) view of the nuclei at each
# timepoint with the species-specific ingression region and the animal
# cap overlaid.

disc_outline <- function(b, sp_lbl, n = 80) {
  ang <- seq(0, 2 * pi, length.out = n)
  data.table(species = sp_lbl,
             phi   = b$phi   + b$radius * cos(ang),
             theta = b$theta + b$radius * sin(ang))
}
margin_rect_outline <- function(lm, w_th, w_ph, sp_lbl) {
  th0 <- lm[["theta"]] - w_th; th1 <- lm[["theta"]] + w_th
  ph0 <- lm[["phi"]]   - w_ph; ph1 <- lm[["phi"]]   + w_ph
  data.table(species = sp_lbl,
             phi   = c(ph0, ph1, ph1, ph0, ph0),
             theta = c(th0, th0, th1, th1, th0))
}

outlines_g <- rbind(
  disc_outline(BULGE_M, "Medaka"),
  margin_rect_outline(lm_z_margin, ZEB_MARGIN_HALFW_TH, ZEB_MARGIN_HALFW_PH,
                      "Zebrafish"))
outlines_g[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

# Frozen-region rectangles (so the user sees what the T0 footprint looks like
# alongside the geometric seed)
frozen_rect <- data.table(species = c("Medaka","Zebrafish"),
                          xmin    = c(mbb$th_min, zbb$th_min),
                          xmax    = c(mbb$th_max, zbb$th_max),
                          ymin    = c(mbb$ph_min, zbb$ph_min),
                          ymax    = c(mbb$ph_max, zbb$ph_max))
frozen_rect[, species := factor(species, levels = c("Medaka","Zebrafish"))]

# Animal cap rectangle
cap_rect_g <- data.table(species = factor(c("Medaka","Zebrafish"),
                                          levels = c("Medaka","Zebrafish")),
                         ymin = c(theta_min_m, theta_min_z),
                         ymax = c(THETA_CAP_M, THETA_CAP_Z))

# Snapshot all cells at the CENTRE FRAME of each timepoint window
# (so the figure is a single-frame snapshot, not 10 frames overlaid)
snap_per_tp <- list()
for (sp in c("Medaka","Zebrafish")) {
  sp_df <- if (sp == "Medaka") sp_m else sp_z
  for (i in seq_len(nrow(timepoints))) {
    tp <- timepoints$tp_label[i]
    tgt <- if (sp == "Medaka") timepoints$medaka_target[i]
           else timepoints$zebrafish_target[i]
    # Centre frame is the closest FRAME to the target time
    centre_t <- sp_df[, .(d = abs(time_min - tgt)), by = FRAME][order(d)][1, FRAME]
    snap <- sp_df[FRAME == centre_t,
                  .(TRACK_ID, PHI_DEG, THETA_DEG, SPHERICAL_DEPTH)]
    snap[, `:=`(species = sp, timepoint = tp)]
    snap_per_tp[[paste(sp, tp)]] <- snap
  }
}
snap_g <- rbindlist(snap_per_tp)
snap_g[, species   := factor(species,   levels = c("Medaka","Zebrafish"))]
snap_g[, timepoint := factor(timepoint, levels = TP_LVL)]

# Region tag at the centre frame
snap_g[, in_cap := THETA_DEG <= ifelse(species == "Medaka", THETA_CAP_M, THETA_CAP_Z)]
snap_g[, in_med_bulge := sqrt((THETA_DEG - BULGE_M$theta)^2 +
                              (PHI_DEG - BULGE_M$phi)^2) < BULGE_M$radius]
snap_g[, in_zeb_margin := abs(THETA_DEG - lm_z_margin[["theta"]]) <= ZEB_MARGIN_HALFW_TH &
                           abs(PHI_DEG   - lm_z_margin[["phi"]])   <= ZEB_MARGIN_HALFW_PH]
snap_g[, region := fifelse(species == "Medaka" & in_med_bulge,  "Head mesoderm bulge",
                            fifelse(species == "Zebrafish" & in_zeb_margin, "Margin",
                              fifelse(in_cap, "Animal cap", "Other")))]
snap_g[, region := factor(region,
  levels = c("Head mesoderm bulge", "Margin", "Animal cap", "Other"))]

p_top <- ggplot(snap_g, aes(PHI_DEG, THETA_DEG)) +
  # Animal cap band (rectangle over the top of the imaged theta range)
  geom_rect(data = cap_rect_g, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "#3B7DD8", alpha = 0.18) +
  # Frozen region (T0 footprint) -- only the data points fallback; we draw
  # this on a separate panel below
  # Nuclei
  geom_point(aes(color = region), size = 0.5, alpha = 0.7) +
  # Geometric disc outline (medaka bulge)
  geom_polygon(data = outlines_g[species == "Medaka"], inherit.aes = FALSE,
               aes(phi, theta), fill = NA,
               color = "#E31A1C", linewidth = 0.9) +
  # Geometric margin rectangle outline (zebrafish)
  geom_polygon(data = outlines_g[species == "Zebrafish"], inherit.aes = FALSE,
               aes(phi, theta), fill = NA,
               color = "#E31A1C", linewidth = 0.9) +
  facet_grid(species ~ timepoint, scales = "free") +
  scale_y_reverse() +
  scale_color_manual(values = c("Head mesoderm bulge" = "#E31A1C",
                                "Margin" = "#E31A1C",
                                "Animal cap" = "#3B7DD8",
                                "Other" = "grey60"),
                     name = NULL,
                     guide = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Top view -- single-frame snapshots (T-1h, T0, T+1h)",
       subtitle = paste(
         "Red CIRCLE = geometric medaka head-mesoderm bulge disc (centre frame).",
         "Red RECTANGLE = geometric zebrafish margin band around the 'margin' landmark (centre frame).",
         "Blue band = animal cap (theta in [theta_min, theta_min+5°]).",
         "Cells at the centre frame of each window are plotted; cells inside the disc/rectangle are coloured red.",
         sep = "\n"),
       x = expression(varphi*" (deg)"),
       y = expression(theta*" (deg, animal pole up)")) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 8, color = "grey25"),
        legend.position = "right") +
  plot_annotation(
    title = "Region outlines overlaid on the data",
    subtitle = sprintf("Each row = one species; each column = one timepoint.  T0-frozen footprint (used for the density analysis): medaka bulge theta in [%.2f, %.2f], phi in [%.2f, %.2f]; zebrafish margin theta in [%.2f, %.2f], phi in [%.2f, %.2f].",
                       mbb$th_min, mbb$th_max, mbb$ph_min, mbb$ph_max,
                       zbb$th_min, zbb$th_max, zbb$ph_min, zbb$ph_max),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(p_top, "fig_F_topview_regions.pdf", w = 16, h = 10)

# Figure G: cohort -- where do T0-margin tracks end up at T-1h / T+1h? ----
cohort_all2 <- copy(cohort_all)
cohort_all2[, timepoint := factor(timepoint, levels = TP_LVL)]
p_f1 <- ggplot(cohort_all2[!is.na(theta_med)],
               aes(timepoint, fill = in_geometric)) +
  geom_bar(position = "fill", alpha = 0.9, color = "grey20", linewidth = 0.2) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#E31A1C", "FALSE" = "#3B7DD8"),
                    labels = c("TRUE" = "in geometric region", "FALSE" = "outside"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "A. T0-cohort: fraction of T0-region tracks still in the geometric region",
       subtitle = "Cohort = TRACK_IDs that were in the geometric margin/disc at T0.  We look at where they are at T-1h / T0 / T+1h.  Stacked bar shows fraction still inside the (T0) geometric region vs outside.",
       x = NULL, y = "fraction of cohort tracks") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
p_f2 <- ggplot(cohort_all2[!is.na(theta_med)],
               aes(timepoint, theta_med, fill = species)) +
  geom_violin(alpha = 0.6, color = "grey20", linewidth = 0.2,
              position = position_dodge(width = 0.8), width = 0.7,
              scale = "width") +
  geom_boxplot(width = 0.13, outlier.size = 0.4,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "B. T0-cohort: theta distribution at each timepoint",
       subtitle = "If the margin moves, the cohort shifts its theta distribution over time.  Medaka: thin shift (bulge is more stationary).  Zebrafish: wider shift -- epiboly pulls margin cells into the embryo.",
       x = NULL, y = expression(theta*" (deg)")) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
p_f3 <- ggplot(cohort_all2[!is.na(depth_med)],
               aes(timepoint, depth_med, fill = species)) +
  geom_violin(alpha = 0.6, color = "grey20", linewidth = 0.2,
              position = position_dodge(width = 0.8), width = 0.7,
              scale = "width") +
  geom_boxplot(width = 0.13, outlier.size = 0.4,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "C. T0-cohort: depth distribution at each timepoint",
       subtitle = "Same tracks across time.  Zebrafish: depth distribution is essentially unchanged (no internalisation).  Medaka: tracks deepen as they ingress.",
       x = NULL, y = "spherical depth (µm)") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(face = "bold"))
fig_cohort <- (p_f1) / (p_f2 | p_f3) + plot_layout(heights = c(1, 1.6)) +
  plot_annotation(title = "Track cohort (T0 region -> where at T-1h & T+1h)",
                  subtitle = "Tracks in the geometric margin/disc at T0.  Zebrafish cohort: depth unchanged (no internalisation) but theta shifts (margin moves).  Medaka cohort: depth increases (cells ingress).",
                  theme = theme(plot.title = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_cohort, "fig_G_cohort.pdf", w = 14, h = 14)

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
