# =============================================================================
# NUCLEI DENSITY @ 3 TIMEPOINTS -- MEDAKA vs ZEBRAFISH
# -----------------------------------------------------------------------------
# ONE-SHOT analysis: only the 3 timepoints T-1h, T0 (ingression), T+1h.
# No continuum time series.
#
# Per timepoint x species x region, computes:
#   * column-integrated density  (n / sphere-surface area, all depths stacked)
#   * volumetric density        (n / area * observed depth range)
#   * KNN distance              (k = 1, 10)
#   * mean NN distance          (= KNN k=1)
#
# Background normalisation: animal cap (theta <= margin - W).
# Non-normalised densities are also shown alongside the fold change.
#
# Time alignment: each species aligned to its own ingression time
#   medaka   : frame 199 = 99.5 min  (FI = 30 s)
#   zebrafish: frame  40 = 80.0 min  (FI = 120 s)
# So "T-1h" = 60 min before each species' own ingression, etc.
#
# Outputs (in results/density_3timepoints/):
#   metrics_long.csv              full long table
#   metrics_wide_with_fold.csv    wide table with fold changes
#   density_3tp_main.pdf          density (per um² & per um³) + fold changes
#   density_3tp_neighbours.pdf    mean NN + KNN k=10
#   density_3tp_validation.pdf    top view of cells + independent density check
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

OUT_DIR <- "results/density_3timepoints"
dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Parameters (copied verbatim from oneshot_comparison.R)
# -----------------------------------------------------------------------------
MEDAKA_INPUT     <- "data/oriented_medaka_ultrack"
ZEBRAFISH_INPUT  <- "data/oriented_zebrafish_ultrack"

MEDAKA_FI  <- 30;  ZEB_FI  <- 120
MEDAKA_VOXEL_UM <- 1.05152;  ZEB_VOXEL_UM <- 1.24785
MEDAKA_INGRESSION_FRAME <- 199L;  ZEB_INGRESSION_FRAME <- 40L

ZONE_MARG_W_M  <- 5    # animal-cap margin half-width (deg)
ZONE_MARG_W_Z  <- 8

# Animal cap is the uppermost (smallest theta) part of the imaged region --
# a "first-latitude" band of width CAP_WIDTH_DEG starting at theta_min.
# We do NOT extend it all the way to the margin; the cap is just a
# representative background sample of the upper cap.
CAP_WIDTH_DEG  <- 5    # width of the upper-cap band (deg)

MEDAKA_MIN_FRAMES <- 20L
ZEB_MIN_FRAMES    <- 5L

INGR_DEPTH_PCTILE <- 0.95
TP_WINDOW_MIN     <- 5     # ±5 min around each target timepoint
KNN_SUBSAMPLE_N   <- 20000 # safety cap for kd-tree queries
KNN_KS            <- c(1L, 10L)

species_colors <- c("Medaka" = "#E69F00", "Zebrafish" = "#0072B2")
REGION_COL <- c("Head mesoderm disc" = "#E31A1C",
                "Animal cap"         = "#3B7DD8")
REGION_LVL <- c("Head mesoderm disc", "Animal cap")
TP_LVL     <- c("T-1h", "T0", "T+1h")

theme_pub <- function(bs = 11) {
  theme_minimal(base_size = bs) +
    theme(plot.title    = element_text(face = "bold", size = bs + 1),
          plot.subtitle = element_text(size = bs - 1, color = "grey40"),
          strip.text    = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey96", color = NA),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom")
}

save_pdf <- function(p, name, w = 12, h = 8) {
  path <- file.path(OUT_DIR, name)
  tmp <- tempfile(pattern = "p_", fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  ggsave(tmp, p, width = w, height = h, device = "pdf")
  ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) stop("failed to write: ", path)
  cat(sprintf("  saved %s\n", path))
}

banner <- function(x) {
  cat("\n", strrep("=", 70), "\n  ", x, "\n", strrep("=", 70), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
banner("LOAD DATA")
sp_m <- fread(file.path(MEDAKA_INPUT, "oriented_tracks_medaka.csv"))
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

cat(sprintf("  Medaka:    %s spots, %d tracks, frames %d-%d\n",
            format(nrow(sp_m), big.mark = ","), uniqueN(sp_m$TRACK_ID),
            min(sp_m$FRAME), max(sp_m$FRAME)))
cat(sprintf("  Zebrafish: %s spots, %d tracks, frames %d-%d\n",
            format(nrow(sp_z), big.mark = ","), uniqueN(sp_z$TRACK_ID),
            min(sp_z$FRAME), max(sp_z$FRAME)))

# -----------------------------------------------------------------------------
# Sphere params, landmarks, data-driven bulge
# -----------------------------------------------------------------------------
read_sphere <- function(dir, vox)
  as.numeric(fread(file.path(dir, "sphere_params.csv"))[parameter == "radius", value]) * vox
R_M <- read_sphere(MEDAKA_INPUT, MEDAKA_VOXEL_UM)
R_Z <- read_sphere(ZEBRAFISH_INPUT, ZEB_VOXEL_UM)

read_landmark <- function(dir, name) {
  d <- fread(file.path(dir, "gastrulation_landmarks.csv"))
  c(theta = as.numeric(d[landmark == name, theta]),
    phi   = as.numeric(d[landmark == name, phi]))
}
MARGIN_M <- read_landmark(MEDAKA_INPUT,    "margin")[["theta"]]
MARGIN_Z <- read_landmark(ZEBRAFISH_INPUT, "margin")[["theta"]]
lm_m_ingr <- read_landmark(MEDAKA_INPUT,    "ingression_center")
lm_z_ingr <- read_landmark(ZEBRAFISH_INPUT, "ingression_center")

detect_bulge <- function(sp, f_theta, f_phi) {
  late <- sp[FRAME >= max(FRAME)/2 & is.finite(SPHERICAL_DEPTH)]
  thr  <- as.numeric(quantile(late$SPHERICAL_DEPTH, INGR_DEPTH_PCTILE, na.rm = TRUE))
  cand <- late[SPHERICAL_DEPTH >= thr]
  if (nrow(cand) >= 50) {
    ct <- mean(cand$THETA_DEG); cp <- mean(cand$PHI_DEG)
    cand[, ang_dist := sqrt((THETA_DEG - ct)^2 + (PHI_DEG - cp)^2)]
    list(theta = ct, phi = cp,
         radius = min(as.numeric(quantile(cand$ang_dist, 0.75, na.rm = TRUE)), 7),
         n_cand = nrow(cand), depth_thresh = thr, source = "data-driven")
  } else {
    list(theta = f_theta, phi = f_phi, radius = 10, n_cand = nrow(cand),
         depth_thresh = thr, source = "landmark fallback")
  }
}
BULGE_M <- detect_bulge(sp_m, lm_m_ingr[["theta"]], lm_m_ingr[["phi"]])
BULGE_Z <- detect_bulge(sp_z, lm_z_ingr[["theta"]], lm_z_ingr[["phi"]])
cat(sprintf("  Bulge: medaka theta=%.1f phi=%.1f r=%.1f | zebrafish theta=%.1f phi=%.1f r=%.1f\n",
            BULGE_M$theta, BULGE_M$phi, BULGE_M$radius,
            BULGE_Z$theta, BULGE_Z$phi, BULGE_Z$radius))

# -----------------------------------------------------------------------------
# Region areas
# -----------------------------------------------------------------------------
# Disc area: a disc in (theta, phi) space of angular radius b$radius
# (degrees) projects onto the sphere with area
#   A = R^2 * sin(theta_center) * pi * (alpha_rad)^2
# where alpha_rad = b$radius * pi/180.  This is the Jacobian-weighted
# small-circle approximation.  (Equivalent to the exact spherical-cap
# formula 2*pi*R^2*(1-cos(alpha)) to leading order in alpha, with the
# sin(theta) factor accounting for the (theta, phi) -> surface area
# Jacobian.)  Kept as-is -- the disc is defined in (theta, phi) space by
# `detect_bulge()`, so this is the consistent area.
disc_area <- function(R, b)
  R^2 * sin(b$theta * pi / 180) * pi * (b$radius * pi / 180)^2
A_BULGE_M <- disc_area(R_M, BULGE_M)
A_BULGE_Z <- disc_area(R_Z, BULGE_Z)

# Animal cap area: integration R^2 sin(theta) dtheta dphi over the IMAGED
# patch (theta from imaged min to tcap, phi from imaged range).  Using the
# full spherical cap (theta from 0 to tcap) over-estimates the area by 2-3x
# because the imaged data only covers theta from ~62 deg (medaka) / ~64 deg
# (zebrafish) up to tcap.  Disc area is computed via the same convention
# but does not need this correction because the disc is small enough to fit
# inside the imaged theta range.
imaged_animal_cap_area <- function(R, theta_cap_deg, theta_min_deg, phi_rng_rad) {
  th_cap <- min(theta_cap_deg, 180) * pi / 180
  th_min <- theta_min_deg * pi / 180
  R^2 * (cos(th_min) - cos(th_cap)) * phi_rng_rad
}
phi_range_rad <- function(dt) {
  q <- as.numeric(quantile(dt$PHI_DEG, c(0.02, 0.98), na.rm = TRUE)) * pi / 180
  q[2] - q[1]
}
theta_min_deg <- function(dt) {
  as.numeric(quantile(dt$THETA_DEG, 0.02, na.rm = TRUE))
}
phi_rng_m <- phi_range_rad(sp_m)
phi_rng_z <- phi_range_rad(sp_z)
theta_min_m <- theta_min_deg(sp_m)
theta_min_z <- theta_min_deg(sp_z)
# Animal cap = upper cap = first-latitude band, theta in [theta_min, theta_min + CAP_WIDTH_DEG]
THETA_CAP_M <- theta_min_m + CAP_WIDTH_DEG
THETA_CAP_Z <- theta_min_z + CAP_WIDTH_DEG
A_ANIM_M <- imaged_animal_cap_area(R_M, THETA_CAP_M, theta_min_m, phi_rng_m)
A_ANIM_Z <- imaged_animal_cap_area(R_Z, THETA_CAP_Z, theta_min_z, phi_rng_z)
cat(sprintf("  Imaged theta min: medaka %.2f, zebrafish %.2f deg\n",
            theta_min_m, theta_min_z))
cat(sprintf("  Imaged phi range: medaka %.2f rad, zebrafish %.2f rad\n",
            phi_rng_m, phi_rng_z))
cat(sprintf("  Areas (um^2, imaged): medaka disc=%.0f, animal cap=%.0f | zebrafish disc=%.0f, animal cap=%.0f\n",
            A_BULGE_M, A_ANIM_M, A_BULGE_Z, A_ANIM_Z))

# -----------------------------------------------------------------------------
# Timepoints
# -----------------------------------------------------------------------------
t_ingr_m <- MEDAKA_INGRESSION_FRAME * MEDAKA_FI / 60
t_ingr_z <- ZEB_INGRESSION_FRAME    * ZEB_FI    / 60
cat(sprintf("  Ingression time: medaka %.1f min, zebrafish %.1f min\n",
            t_ingr_m, t_ingr_z))

timepoints <- data.table(
  tp_label = TP_LVL,
  offset_min = c(-60, 0, 60))
timepoints[, medaka_target    := t_ingr_m + offset_min]
timepoints[, zebrafish_target := t_ingr_z + offset_min]
print(timepoints)

# -----------------------------------------------------------------------------
# Selectors + per-region metrics
# -----------------------------------------------------------------------------
in_disc        <- function(sp, b) sp[!is.na(PHI_DEG) &
  sqrt((THETA_DEG - b$theta)^2 + (PHI_DEG - b$phi)^2) < b$radius]
in_animal_cap  <- function(sp, theta_cap) sp[!is.na(THETA_DEG) & THETA_DEG <= theta_cap]

knn_means <- function(coords, ks = KNN_KS, cap = KNN_SUBSAMPLE_N) {
  n <- nrow(coords)
  if (n <= max(ks)) return(setNames(rep(NA_real_, length(ks)), paste0("k", ks)))
  if (n > cap) {
    set.seed(1L)
    coords <- coords[sample.int(n, cap), , drop = FALSE]
  }
  nn <- RANN::nn2(coords, k = max(ks) + 1L)
  d  <- nn$nn.dists[, -1, drop = FALSE]
  setNames(vapply(ks, function(k)
                    if (k <= ncol(d)) mean(d[, k]) else NA_real_, numeric(1)),
           paste0("k", ks))
}

region_metrics <- function(cells, area_um2) {
  if (nrow(cells) == 0)
    return(list(n_cells = 0L, area_um2 = area_um2,
                depth_range_um = NA_real_,
                volume_um3 = NA_real_, density_per_um2 = NA_real_,
                density_per_um3 = NA_real_, knn_k1_um = NA_real_,
                knn_k10_um = NA_real_, mean_nn_um = NA_real_))
  # Use ONE position per track (median over the time window) so we don't
  # double-count cells that appear at multiple timepoints within the window.
  # Without this, n_cells is ~14x too high and the KNN is ~14x too small
  # (the nearest "neighbour" turns out to be the same cell at the previous
  # frame, ~0.5-1 um away).
  cells_uniq <- cells[, .(
    POSITION_X     = median(POSITION_X,     na.rm = TRUE),
    POSITION_Y     = median(POSITION_Y,     na.rm = TRUE),
    POSITION_Z     = median(POSITION_Z,     na.rm = TRUE),
    SPHERICAL_DEPTH = median(SPHERICAL_DEPTH, na.rm = TRUE)
  ), by = TRACK_ID]
  n <- nrow(cells_uniq)
  if (n == 0)
    return(list(n_cells = 0L, area_um2 = area_um2,
                depth_range_um = NA_real_,
                volume_um3 = NA_real_, density_per_um2 = NA_real_,
                density_per_um3 = NA_real_, knn_k1_um = NA_real_,
                knn_k10_um = NA_real_, mean_nn_um = NA_real_))
  # 2nd-98th percentile depth range (robust to outliers in SPHERICAL_DEPTH)
  depth <- cells_uniq[is.finite(SPHERICAL_DEPTH), SPHERICAL_DEPTH]
  depth_range <- if (length(depth) >= 2)
    as.numeric(quantile(depth, 0.98) - quantile(depth, 0.02)) else 0
  volume_um3 <- area_um2 * depth_range
  coords <- as.matrix(cells_uniq[, .(POSITION_X, POSITION_Y, POSITION_Z)])
  knn <- knn_means(coords)
  list(
    n_cells        = n,
    area_um2       = area_um2,
    depth_range_um = depth_range,
    volume_um3     = volume_um3,
    density_per_um2 = n / area_um2,
    density_per_um3 = if (volume_um3 > 0) n / volume_um3 else NA_real_,
    knn_k1_um      = unname(knn["k1"]),
    knn_k10_um     = unname(knn["k10"]),
    mean_nn_um     = unname(knn["k1"]))
}

# -----------------------------------------------------------------------------
# Compute metrics at the 3 timepoints
# -----------------------------------------------------------------------------
banner("COMPUTE METRICS @ 3 TIMEPOINTS")

sp_list   <- list(Medaka = sp_m, Zebrafish = sp_z)
bulge_l   <- list(Medaka = BULGE_M, Zebrafish = BULGE_Z)
A_disc_l  <- list(Medaka = A_BULGE_M, Zebrafish = A_BULGE_Z)
A_anim_l  <- list(Medaka = A_ANIM_M, Zebrafish = A_ANIM_Z)
tcap_l    <- list(Medaka = THETA_CAP_M, Zebrafish = THETA_CAP_Z)
target_l  <- list(Medaka = timepoints$medaka_target,
                  Zebrafish = timepoints$zebrafish_target)

rows <- list()
for (sp_lbl in c("Medaka", "Zebrafish")) {
  sp     <- sp_list[[sp_lbl]]
  b      <- bulge_l[[sp_lbl]]
  for (i in seq_len(nrow(timepoints))) {
    tgt <- target_l[[sp_lbl]][i]
    win <- sp[time_min >= tgt - TP_WINDOW_MIN & time_min <= tgt + TP_WINDOW_MIN]
    rows[[length(rows) + 1L]] <- cbind(
      data.table(species = sp_lbl, timepoint = timepoints$tp_label[i],
                 offset_min = timepoints$offset_min[i],
                 target_time_min = tgt, region = "Head mesoderm disc"),
      as.data.table(region_metrics(in_disc(win, b), A_disc_l[[sp_lbl]])))
    rows[[length(rows) + 1L]] <- cbind(
      data.table(species = sp_lbl, timepoint = timepoints$tp_label[i],
                 offset_min = timepoints$offset_min[i],
                 target_time_min = tgt, region = "Animal cap"),
      as.data.table(region_metrics(in_animal_cap(win, tcap_l[[sp_lbl]]),
                                    A_anim_l[[sp_lbl]])))
  }
}
metrics <- rbindlist(rows)
metrics[, species   := factor(species,   levels = c("Medaka", "Zebrafish"))]
metrics[, region    := factor(region,    levels = REGION_LVL)]
metrics[, timepoint := factor(timepoint, levels = TP_LVL)]

# Add context: expected NN if cells were uniformly distributed in 3D
# (matching the KNN measurement, which is in 3D).
#   3D Poisson expected NN = (3 / (4*pi*density_per_um3))^(1/3)
#   clustering_ratio = actual / expected.
#   ratio < 1 -> cells are clustered (locally denser than uniform 3D Poisson).
#   ratio ~ 1 -> uniform.  ratio > 1 -> hyper-dispersed.
metrics[, expected_nn_3D_uniform_um :=
          (3 / (4 * pi * density_per_um3))^(1/3)]
metrics[, clustering_ratio := mean_nn_um / expected_nn_3D_uniform_um]

# Wide form with fold change (head mesoderm / animal cap)
metrics_wide <- dcast(
  metrics, species + timepoint + offset_min + target_time_min ~ region,
  value.var = c("n_cells",
                "density_per_um2", "density_per_um3",
                "knn_k1_um", "knn_k10_um", "mean_nn_um",
                "expected_nn_3D_uniform_um", "clustering_ratio"))
metrics_wide[, fold_density_um2 := `density_per_um2_Head mesoderm disc` /
                                      `density_per_um2_Animal cap`]
metrics_wide[, fold_density_um3 := `density_per_um3_Head mesoderm disc` /
                                      `density_per_um3_Animal cap`]

fwrite(metrics,      file.path(OUT_DIR, "metrics_long.csv"))
fwrite(metrics_wide, file.path(OUT_DIR, "metrics_wide_with_fold.csv"))
cat("\n  Console summary (clustering ratio = mean_nn / 3D-Poisson-expected_nn):\n")
print(metrics[, .(species, timepoint, region, n_cells,
                  depth_um = round(depth_range_um, 1),
                  density_per_um2 = round(density_per_um2, 4),
                  density_per_um3 = round(density_per_um3, 5),
                  nn_um = round(mean_nn_um, 2),
                  expected_nn_3D = round(expected_nn_3D_uniform_um, 2),
                  clustering_ratio = round(clustering_ratio, 2))])

# -----------------------------------------------------------------------------
# FIGURE 1: NUCLEI DENSITY  --  real units (per um^2, per um^3) + disc/cap fold
# -----------------------------------------------------------------------------
banner("FIGURE 1 -- nuclei density")

# Bar plot of a density metric, faceted by region.  Below each bar we
# annotate the sample (n_cells) and -- for volumetric panels -- the actual
# volume the cells were counted in.
bar_density <- function(metrics, ycol, ylab, title, subtitle,
                        ann2_col = NULL, ann2_fmt = NULL) {
  d <- copy(metrics)
  d[, y := get(ycol)]
  p <- ggplot(d, aes(timepoint, y, fill = species)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7,
             alpha = 0.9, color = "grey20", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.3g", y)),
              position = position_dodge(width = 0.78),
              vjust = -0.4, size = 3) +
    geom_text(aes(label = paste0("n=", format(n_cells, big.mark = ","))),
              position = position_dodge(width = 0.78),
              vjust = 1.3, size = 2.2, color = "grey25")
  if (!is.null(ann2_col)) {
    p <- p + geom_text(aes(label = sprintf(ann2_fmt, get(ann2_col))),
                       position = position_dodge(width = 0.78),
                       vjust = 2.8, size = 2.4, color = "grey30")
  }
  p + facet_wrap(~ region, scales = "free_y") +
    scale_fill_manual(values = species_colors, name = NULL) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 8, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
}

# Bar plot of a fold change (disc / cap), single value per species/timepoint.
# 1.0 line for "no enrichment".
bar_fold <- function(metrics_wide, ycol, ylab, title, subtitle) {
  d <- metrics_wide[, .(species, timepoint, y = get(ycol))]
  d[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
  ggplot(d, aes(timepoint, y, fill = species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             alpha = 0.9, color = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
    geom_text(aes(label = sprintf("%.2f", y)),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3.4) +
    scale_fill_manual(values = species_colors, name = NULL) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
}

# Panel A: density per µm².  Surface = full disc area (a 2D circle in
# (theta, phi)) or full cap area (a 5° band at the top).  All nuclei in
# the region are divided by this surface.
p_d_um2 <- bar_density(
  metrics, "density_per_um2",
  expression("nuclei / "*mu*"m"^2),
  "A. Density per µm²",
  sprintf(
    "Full surface area (µm², constant per region):  medaka disc %s / cap %s ;  zebrafish disc %s / cap %s.  n / surface.",
    format(round(A_BULGE_M), big.mark = ","), format(round(A_ANIM_M), big.mark = ","),
    format(round(A_BULGE_Z), big.mark = ","), format(round(A_ANIM_Z), big.mark = ",")))

# Panel B: density per µm³.  Volume = full surface × full depth range.
# Shown under each bar in µm³ (full number, comma-separated).
metrics[, volume_um3_str := format(round(volume_um3), big.mark = ",")]
p_d_um3 <- bar_density(
  metrics, "density_per_um3",
  expression("nuclei / "*mu*"m"^3),
  "B. Density per µm³",
  "Volume = full surface × full depth range.  vol shown below each bar in µm³.",
  ann2_col = "volume_um3_str", ann2_fmt = "vol=%s µm³")

# Panel C: fold per µm² (disc / cap)
p_fold_um2 <- bar_fold(
  metrics_wide, "fold_density_um2",
  "disc / cap",
  "C. Fold per µm²  (disc / cap)",
  "Disc density divided by cap density.  1.0 = no enrichment.")

# Panel D: fold per µm³ (disc / cap)
p_fold_um3 <- bar_fold(
  metrics_wide, "fold_density_um3",
  "disc / cap",
  "D. Fold per µm³  (disc / cap)",
  "Same as C, volumetric (comparable across species).")

fig_density <- (p_d_um2 | p_d_um3) / (p_fold_um2 | p_fold_um3) +
  plot_annotation(
    title = "Nuclei density  --  Head mesoderm disc vs Animal cap (5° band)",
    subtitle = sprintf(
      "T0 = species-specific ingression (medaka frame %d = %.1f min; zebrafish frame %d = %.1f min).  Per-µm² stacks every z-depth; per-µm³ divides by depth range.  Fold = disc / cap.",
      MEDAKA_INGRESSION_FRAME, t_ingr_m,
      ZEB_INGRESSION_FRAME,    t_ingr_z),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_density, "density_3tp_metrics.pdf", w = 16, h = 11)

# -----------------------------------------------------------------------------
# FIGURE 2: NUCLEI COUNTS + LOCAL PACKING  --  what the density numbers rest on
# -----------------------------------------------------------------------------
banner("FIGURE 2 -- counts and local packing")

# Panel A: n_cells (line) per region over time.  The numerator of every
# density.  Lines connect the same species across timepoints.
plot_n_cells <- function(metrics) {
  ggplot(metrics, aes(timepoint, n_cells, color = species, group = species)) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 3) +
    geom_text(aes(label = format(n_cells, big.mark = ",")),
              vjust = -0.8, size = 2.8, show.legend = FALSE) +
    facet_wrap(~ region, scales = "free_y") +
    scale_color_manual(values = species_colors, name = NULL) +
    labs(title = "A. n_cells per region",
         subtitle = "Numerator of every density.  Medaka cap shrinks; medaka disc and zebrafish both grow.",
         x = NULL, y = "nuclei in window") +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
}

# Panel B: depth_range (tissue column thickness) per region over time.
# The extra factor for per-µm³ (vs per-µm²).
plot_depth_range <- function(metrics) {
  ggplot(metrics, aes(timepoint, depth_range_um, fill = species)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7,
             alpha = 0.9, color = "grey20", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.0f", depth_range_um)),
              position = position_dodge(width = 0.78),
              vjust = -0.4, size = 3) +
    facet_wrap(~ region, scales = "free_y") +
    scale_fill_manual(values = species_colors, name = NULL) +
    labs(title = "B. Depth range  (µm)",
         subtitle = "Tissue column thickness.  Deeper columns mean more vertical room for nuclei to stack.",
         x = NULL, y = expression("depth range ("*mu*"m)")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
}

# Panels C, D: KNN bar plot WITH 3D-Poisson expected-NN reference line.
# Bar BELOW the dashed line -> cells are CLUSTERED.
# Bar ABOVE the dashed line -> cells are HYPER-DISPERSED.
# Bar AT the dashed line -> uniform (3D Poisson).
# This is more direct than a ratio: it shows real spacing in µm.
knn_bar_with_expected <- function(metrics, ycol, ylab, title, subtitle) {
  d <- copy(metrics)
  d[, y := get(ycol)]
  ggplot(d, aes(timepoint, y, fill = species)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7,
             alpha = 0.9, color = "grey20", linewidth = 0.2) +
    # 3D-Poisson expected NN as a per-bar dashed tick
    geom_segment(aes(x = as.numeric(timepoint) - 0.30,
                     xend = as.numeric(timepoint) + 0.30,
                     y = expected_nn_3D_uniform_um,
                     yend = expected_nn_3D_uniform_um,
                     color = species),
                 position = position_dodge(width = 0.78),
                 linewidth = 0.5, linetype = "dashed", show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f", y)),
              position = position_dodge(width = 0.78),
              vjust = -0.5, size = 2.8) +
    facet_wrap(~ region, scales = "free_y") +
    scale_fill_manual(values = species_colors, name = NULL) +
    scale_color_manual(values = species_colors, name = NULL) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x = element_text(face = "bold"))
}

# Panel C: KNN k=1
p_knn1 <- knn_bar_with_expected(
  metrics, "mean_nn_um",
  expression("3D distance to nearest nucleus ("*mu*"m)"),
  "C. KNN k=1  (3D distance to nearest nucleus)",
  "Real spacing in µm.  Dashed line = expected distance for a uniform 3D distribution at the same volumetric density.")

# Panel D: KNN k=10
p_knn10 <- knn_bar_with_expected(
  metrics, "knn_k10_um",
  expression("3D distance to 10th nucleus ("*mu*"m)"),
  "D. KNN k=10  (3D distance to 10th-closest nucleus)",
  "Same idea, less sensitive to single-cell outliers.")

fig_packing <- (plot_n_cells(metrics) | plot_depth_range(metrics)) /
                (p_knn1            | p_knn10) +
  plot_annotation(
    title = "Counts and local packing  --  what the density numbers rest on",
    subtitle = "A and B are the inputs to density (per µm² = n / area; per µm³ = n / (area * depth)).  C and D show local packing directly in µm, with the 3D-Poisson expectation drawn as a dashed tick for each bar (bar below the tick = clustered; bar above = hyper-dispersed).",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_packing, "density_3tp_drivers.pdf", w = 16, h = 11)

# -----------------------------------------------------------------------------
# FIGURE 3: VALIDATION -- show real data + independent density check
# -----------------------------------------------------------------------------
banner("FIGURE 3 -- validation")

# Region overlay helpers
disc_outline <- function(b, sp_lbl, n = 80) {
  ang <- seq(0, 2 * pi, length.out = n)
  data.table(species = sp_lbl,
             phi   = b$phi   + b$radius * cos(ang),
             theta = b$theta + b$radius * sin(ang))
}
outlines <- rbind(disc_outline(BULGE_M, "Medaka"),
                  disc_outline(BULGE_Z, "Zebrafish"))
outlines[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

# Cap is a spherical cap centred on the animal pole (phi = pole meridian,
# theta = 0).  In the (phi, theta) projection the cap projects to a CIRCLE
# of radius theta_cap centred at (phi_animal_pole, 0).  We only image the
# lower portion of that circle (theta >= theta_min), so we trace a closed
# polygon with a flat top at theta_min and two arcs (the cap boundary) on
# the sides -- a "lens" / circular segment.  This is a true cap, NOT a
# rectangular band spanning the whole embryo's phi range.
cap_lens <- function(theta_cap, theta_min, phi_center, sp_lbl, n = 120) {
  stopifnot(theta_min < theta_cap)
  phi_max <- sqrt(theta_cap^2 - theta_min^2)
  t_start <- asin(theta_min / theta_cap)
  # Right arc: from (phi_center + phi_max, theta_min) down to (phi_center, theta_cap)
  t_right <- seq(t_start, pi / 2, length.out = n)
  right_arc <- data.table(phi = phi_center + theta_cap * cos(t_right),
                          theta =              theta_cap * sin(t_right))
  # Left arc: from (phi_center, theta_cap) up to (phi_center - phi_max, theta_min)
  t_left  <- seq(pi / 2, pi - t_start, length.out = n)
  left_arc  <- data.table(phi = phi_center + theta_cap * cos(t_left),
                          theta =              theta_cap * sin(t_left))
  # Top edge: horizontal at theta_min (two endpoints)
  top <- data.table(phi = c(phi_center + phi_max, phi_center - phi_max),
                    theta = theta_min)
  rbind(top, right_arc, left_arc)[, species := sp_lbl]
}
# Animal cap as a thin horizontal band at the TOP of the plot (the "first
# latitude" of the embryo, theta in [theta_min, theta_min + CAP_WIDTH_DEG]).
# This is the upper cap, NOT a band reaching to the margin.
anim_band_rect <- data.table(
  species = factor(c("Medaka", "Zebrafish"),
                   levels = c("Medaka", "Zebrafish")),
  ymin    = c(theta_min_m, theta_min_z),
  ymax    = c(THETA_CAP_M, THETA_CAP_Z))

# Disc band rectangle (for the 1D theta histogram, Panel B) -- shows the
# disc's theta range on the histogram.
disc_band <- data.table(
  species = factor(c("Medaka", "Zebrafish"),
                   levels = c("Medaka", "Zebrafish")),
  xmin    = c(BULGE_M$theta - BULGE_M$radius, BULGE_Z$theta - BULGE_Z$radius),
  xmax    = c(BULGE_M$theta + BULGE_M$radius, BULGE_Z$theta + BULGE_Z$radius))

BIN <- 4  # deg per bin for validation grid

# Per-frame region tag (for any timepoint)
tag_region <- function(sp, b, tcap) {
  sp[, in_disc := sqrt((THETA_DEG - b$theta)^2 + (PHI_DEG - b$phi)^2) < b$radius]
  sp[, in_anim := THETA_DEG <= tcap]
  sp[, region := fifelse(in_disc, "Head mesoderm disc",
                         fifelse(in_anim, "Animal cap", "Other"))]
  sp
}
sp_m_tagged <- tag_region(sp_m, BULGE_M, THETA_CAP_M)
sp_z_tagged <- tag_region(sp_z, BULGE_Z, THETA_CAP_Z)

# -----------------------------------------------------------------------------
# PANEL A: Cells at each timepoint, side by side (top view)
# -----------------------------------------------------------------------------
build_top_view <- function(sp_tagged, sp_lbl, target_t, tp_lbl) {
  win <- sp_tagged[time_min >= target_t - TP_WINDOW_MIN &
                   time_min <= target_t + TP_WINDOW_MIN]
  win[, species := sp_lbl]
  win[, timepoint := tp_lbl]
  win
}
snap_list <- list(
  build_top_view(sp_m_tagged, "Medaka",    t_ingr_m - 60, "T-1h"),
  build_top_view(sp_m_tagged, "Medaka",    t_ingr_m,      "T0"),
  build_top_view(sp_m_tagged, "Medaka",    t_ingr_m + 60, "T+1h"),
  build_top_view(sp_z_tagged, "Zebrafish", t_ingr_z - 60, "T-1h"),
  build_top_view(sp_z_tagged, "Zebrafish", t_ingr_z,      "T0"),
  build_top_view(sp_z_tagged, "Zebrafish", t_ingr_z + 60, "T+1h"))
snap_all <- rbindlist(snap_list)
snap_all[, species   := factor(species,   levels = c("Medaka", "Zebrafish"))]
snap_all[, timepoint := factor(timepoint, levels = TP_LVL)]
snap_all[, region    := factor(region,    levels = c("Head mesoderm disc",
                                                       "Animal cap", "Other"))]

p_top <- ggplot(snap_all, aes(PHI_DEG, THETA_DEG)) +
  # Animal cap as a thin horizontal band at the TOP of the plot (the
  # "first latitude" of the embryo, theta in [theta_min, theta_min + 5°]).
  # This is the UPPER cap only -- not a band reaching to the margin.
  geom_rect(data = anim_band_rect, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = REGION_COL[["Animal cap"]], alpha = 0.18) +
  geom_polygon(data = outlines, inherit.aes = FALSE,
               aes(phi, theta), fill = NA,
               color = REGION_COL[["Head mesoderm disc"]], linewidth = 0.9) +
  geom_point(aes(color = region), size = 0.45, alpha = 0.7) +
  facet_grid(species ~ timepoint, scales = "free") +
  scale_y_reverse() +
  scale_color_manual(values = c("Head mesoderm disc" = "#E31A1C",
                                "Animal cap"         = "#3B7DD8",
                                "Other"              = "grey60"),
                     name = NULL,
                     guide = guide_legend(override.aes = list(size = 3))) +
  labs(title = "A. Cells at the 3 timepoints (top view) -- where are they, exactly?",
       subtitle = sprintf("Each row = one species. Each column = one timepoint. Blue CIRCLE = the IMAGED part of the animal cap (a circular cap around the animal pole, clipped to theta >= %.1f deg). Red contour = head-mesoderm disc. Cells coloured blue = inside animal cap, red = inside disc, grey = outside both.", theta_min_m),
       x = expression(varphi*" (deg)"),
       y = expression(theta*" (deg, animal pole up)")) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        legend.position = "right")

# -----------------------------------------------------------------------------
# PANEL B: 1D histogram of cells vs theta -- spatial redistribution evidence
# -----------------------------------------------------------------------------
THETA_BIN <- 2
build_theta_hist <- function(sp_tagged, sp_lbl, target_t) {
  win <- sp_tagged[time_min >= target_t - TP_WINDOW_MIN &
                   time_min <= target_t + TP_WINDOW_MIN]
  win[, theta_bin := floor(THETA_DEG / THETA_BIN) * THETA_BIN + THETA_BIN/2]
  win[, .(n = .N), by = .(theta_bin, region)
      ][, .(species = sp_lbl, theta_bin, region, n)]
}
# Each timepoint uses its OWN time window (NOT all three labelled T-1h).
# Previous version had a bug where the first rbindlist used c(t-60, t, t+60)
# but labelled every chunk "T-1h", tripling the T-1h counts.
theta_hist <- rbindlist(c(
  lapply(c(t_ingr_m - 60), function(t) cbind(timepoint = "T-1h", build_theta_hist(sp_m_tagged, "Medaka",    t))),
  lapply(c(t_ingr_z - 60), function(t) cbind(timepoint = "T-1h", build_theta_hist(sp_z_tagged, "Zebrafish", t))),
  lapply(c(t_ingr_m),      function(t) cbind(timepoint = "T0",   build_theta_hist(sp_m_tagged, "Medaka",    t))),
  lapply(c(t_ingr_z),      function(t) cbind(timepoint = "T0",   build_theta_hist(sp_z_tagged, "Zebrafish", t))),
  lapply(c(t_ingr_m + 60), function(t) cbind(timepoint = "T+1h", build_theta_hist(sp_m_tagged, "Medaka",    t))),
  lapply(c(t_ingr_z + 60), function(t) cbind(timepoint = "T+1h", build_theta_hist(sp_z_tagged, "Zebrafish", t)))))
theta_hist[, species   := factor(species,   levels = c("Medaka", "Zebrafish"))]
theta_hist[, timepoint := factor(timepoint, levels = TP_LVL)]
# Stacking order (bottom -> top): Other, Animal cap, Head mesoderm disc.
# This puts the disc (the region of interest) on TOP so its small
# contribution is visible above the cap, not buried at the bottom.
theta_hist[, region    := factor(region,    levels = c("Other", "Animal cap", "Head mesoderm disc"))]

p_hist <- ggplot(theta_hist, aes(theta_bin, n, fill = region)) +
  # Shaded band showing the head-mesoderm disc's theta range on the x axis.
  geom_rect(data = disc_band, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = REGION_COL[["Head mesoderm disc"]], alpha = 0.10) +
  # STACKED bars (default geom_col stacking).  Disc is on top so the red
  # contribution is visible.  Note: disc cells are at the SAME depth range
  # as cap cells (z coverage is complete -- both regions include the full
  # column of nuclei at their (phi, theta) location), they are just a
  # small spatial subset, so they appear as a thin red layer at the top of
  # the cap-only stack in each disc bin.
  geom_col(alpha = 0.9, color = "white", linewidth = 0.15) +
  # Vertical line at animal-cap boundary
  geom_vline(data = data.table(
    species = factor(c("Medaka", "Zebrafish"), levels = c("Medaka", "Zebrafish")),
    xint = c(THETA_CAP_M, THETA_CAP_Z)),
    aes(xintercept = xint), linetype = "dashed", color = "grey30",
    inherit.aes = FALSE) +
  facet_grid(species ~ timepoint, scales = "free_y") +
  scale_fill_manual(values = c("Head mesoderm disc" = "#E31A1C",
                                "Animal cap"         = "#3B7DD8",
                                "Other"              = "grey60"),
                     name = NULL) +
  labs(title = "B. 1D histogram of cells per theta band (2° bins, stacked)",
       subtitle = sprintf("Cells per 2° theta band, stacked by region. Light-red shaded band = head-mesoderm disc's theta range (%.1f-%.1f° medaka, %.1f-%.1f° zebrafish). Dashed line = animal-cap boundary (theta = %.1f° medaka, %.1f° zebrafish). The disc covers the FULL z depth (verified: -12 to 104 µm, same as the cap), so the small red layer is not a depth artefact -- it is just the small spatial extent of the disc circle.",
                         BULGE_M$theta - BULGE_M$radius, BULGE_M$theta + BULGE_M$radius,
                         BULGE_Z$theta - BULGE_Z$radius, BULGE_Z$theta + BULGE_Z$radius,
                         THETA_CAP_M, THETA_CAP_Z),
       x = expression(theta*" (deg, animal pole up)"),
       y = "nuclei in band") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 8, color = "grey25"),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# -----------------------------------------------------------------------------
# PANEL C: TRACK-LEVEL analysis -- how many cells LEFT the animal cap between
# consecutive timepoints?  This is the empirical test of the interpretation.
# A track is in AC at timepoint X if its position at the centre of window X
# has THETA_DEG <= tcap.  Then we count: persistent, left, entered.
# -----------------------------------------------------------------------------
track_in_region <- function(sp_tagged, b, tcap, target_t) {
  # Use the median position of each track WITHIN the window centered at target_t
  win <- sp_tagged[time_min >= target_t - TP_WINDOW_MIN &
                   time_min <= target_t + TP_WINDOW_MIN,
                   .(TRACK_ID, THETA_DEG, PHI_DEG, time_min, region)]
  win[, .(theta_med = median(THETA_DEG, na.rm = TRUE),
          phi_med   = median(PHI_DEG,   na.rm = TRUE)),
      by = TRACK_ID][, in_disc := sqrt((theta_med - b$theta)^2 +
                                        (phi_med - b$phi)^2) < b$radius]
}
track_state_m_T1 <- track_in_region(sp_m_tagged, BULGE_M, THETA_CAP_M, t_ingr_m - 60)
track_state_m_T0 <- track_in_region(sp_m_tagged, BULGE_M, THETA_CAP_M, t_ingr_m)
track_state_m_T2 <- track_in_region(sp_m_tagged, BULGE_M, THETA_CAP_M, t_ingr_m + 60)
track_state_z_T1 <- track_in_region(sp_z_tagged, BULGE_Z, THETA_CAP_Z, t_ingr_z - 60)
track_state_z_T0 <- track_in_region(sp_z_tagged, BULGE_Z, THETA_CAP_Z, t_ingr_z)
track_state_z_T2 <- track_in_region(sp_z_tagged, BULGE_Z, THETA_CAP_Z, t_ingr_z + 60)

# Build transitions: AC_track_ids at T1, AC_track_ids at T0, etc.
ac_m_t1 <- track_state_m_T1[in_disc == FALSE, TRACK_ID]
ac_m_t0 <- track_state_m_T0[in_disc == FALSE, TRACK_ID]
ac_m_t2 <- track_state_m_T2[in_disc == FALSE, TRACK_ID]
ac_z_t1 <- track_state_z_T1[in_disc == FALSE, TRACK_ID]
ac_z_t0 <- track_state_z_T0[in_disc == FALSE, TRACK_ID]
ac_z_t2 <- track_state_z_T2[in_disc == FALSE, TRACK_ID]

transitions <- data.table(
  species   = rep(c("Medaka", "Zebrafish"), each = 4),
  transition = rep(c("AC at T-1h & T0 (persistent)",
                     "AC at T-1h, NOT at T0 (left)",
                     "NOT AC at T-1h, AC at T0 (entered)",
                     "AC at T0 & T+1h (persistent)"), 2),
  n_tracks = c(
    length(intersect(ac_m_t1, ac_m_t0)),
    length(setdiff(ac_m_t1, ac_m_t0)),
    length(setdiff(ac_m_t0, ac_m_t1)),
    length(intersect(ac_m_t0, ac_m_t2)),
    length(intersect(ac_z_t1, ac_z_t0)),
    length(setdiff(ac_z_t1, ac_z_t0)),
    length(setdiff(ac_z_t0, ac_z_t1)),
    length(intersect(ac_z_t0, ac_z_t2))))
transitions[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

# Save
fwrite(transitions, file.path(OUT_DIR, "validation_track_transitions.csv"))
cat("\nTrack transitions in the ANIMAL CAP between consecutive timepoints:\n")
print(transitions)

# Bar plot of transitions
p_trans <- ggplot(transitions,
                  aes(transition, n_tracks, fill = species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           alpha = 0.9, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n_tracks),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = species_colors, name = NULL) +
  labs(title = "C. Track-level transitions in the animal cap",
       subtitle = "Tracks are labelled by their median position WITHIN each timepoint window.  'persistent' = track was in animal cap at both consecutive timepoints.  'left' = was in AC at T-1h but not at T0.  'entered' = was not in AC at T-1h but was at T0.  A large 'left' column is the direct evidence that cells leave the animal cap between timepoints.",
       x = NULL, y = "n tracks") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        axis.text.x = element_text(angle = 12, hjust = 1, size = 9),
        legend.position = "right")

# -----------------------------------------------------------------------------
# Independent binned density check (kept as before but now panel D)
# -----------------------------------------------------------------------------
build_binned_check <- function(sp_tagged, sp_lbl, b, R, target_t) {
  win <- sp_tagged[time_min >= target_t - TP_WINDOW_MIN &
                   time_min <= target_t + TP_WINDOW_MIN]
  if (nrow(win) == 0) return(NULL)
  bin_area <- function(th_deg) R^2 * (BIN * pi/180)^2 * sin(th_deg * pi/180)
  win[, theta_bin := floor(THETA_DEG / BIN) * BIN + BIN/2]
  win[, phi_bin   := floor(PHI_DEG   / BIN) * BIN + BIN/2]
  by_bin <- win[, .(n = .N), by = .(theta_bin, phi_bin)]
  by_bin[, area_um2 := bin_area(theta_bin)]
  by_bin[, density_per_um2 := n / area_um2]
  by_bin[, in_disc := sqrt((theta_bin - b$theta)^2 +
                            (phi_bin   - b$phi)^2)   < b$radius]
  by_bin[, in_anim := theta_bin <= THETA_CAP_M]
  by_bin[, species := sp_lbl]
  by_bin
}
chk_m <- build_binned_check(sp_m_tagged, "Medaka",    BULGE_M, R_M, t_ingr_m)
chk_z <- build_binned_check(sp_z_tagged, "Zebrafish", BULGE_Z, R_Z, t_ingr_z)
chk   <- rbind(chk_m, chk_z)
chk[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

headline <- metrics_wide[timepoint == "T0",
                         .(species, region_disc = `density_per_um2_Head mesoderm disc`,
                           region_anim = `density_per_um2_Animal cap`)]
binned_summary <- chk[, .(
  binned_disc_density = sum(n[in_disc])   / sum(area_um2[in_disc]),
  binned_anim_density = sum(n[in_anim])   / sum(area_um2[in_anim]),
  n_bins_disc = sum(in_disc),
  n_bins_anim = sum(in_anim)
), by = species]
binned_summary <- merge(binned_summary, headline, by = "species")
binned_summary[, ratio_disc := binned_disc_density / region_disc]
binned_summary[, ratio_anim := binned_anim_density / region_anim]
cat("\nIndependent binned-density check at T0 (should match headline ~1.0):\n")
print(binned_summary)
fwrite(binned_summary, file.path(OUT_DIR, "validation_binned_check_T0.csv"))

p_check <- ggplot(chk, aes(phi_bin, theta_bin, fill = density_per_um2)) +
  geom_tile() +
  geom_polygon(data = outlines, inherit.aes = FALSE,
               aes(phi, theta), fill = NA,
               color = "white", linewidth = 0.5) +
  geom_rect(data = anim_band_rect, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = NA, color = "white", linewidth = 0.4,
            linetype = "dashed") +
  facet_wrap(~ species, scales = "free") +
  scale_fill_viridis_c(option = "magma",
                       name = expression("cells / "*mu*"m"^2),
                       trans = "sqrt") +
  scale_y_reverse() +
  labs(title = "D. Independent density check (4° × 4° bins) at T0",
       subtitle = sprintf("Per-bin density (cells / µm²) on a %d° x %d° grid. White dashed = animal-cap circle (imaged part). Should match headline density.", BIN, BIN),
       x = expression(varphi*" (deg)"),
       y = expression(theta*" (deg)")) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"))

# Save the validation PDF (page 1: cells at each timepoint + transitions)
fig_val1 <- (p_top / p_hist / p_trans) +
  plot_annotation(
    title = "VALIDATION -- where are the cells, do they really leave the animal cap?",
    subtitle = "A: top view per timepoint, blue band is the IMAGED animal cap (NOT -Inf). B: 1D theta histogram, blue bars should shrink at high theta between timepoints if cells leave. C: track-level transitions -- 'left' counts the tracks that were in AC at T-1h and are no longer at T0. This is the empirical test of the interpretation.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, color = "grey25")))
save_pdf(fig_val1, "density_3tp_validation.pdf", w = 16, h = 22)

# Save the binned check as a separate small PDF (page 2 of validation)
save_pdf(p_check, "density_3tp_validation_binned.pdf", w = 16, h = 8)

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))