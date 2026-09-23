# =============================================================================
# Head mesoderm density normalised by the ANIMAL POLE (animal cap)
# -----------------------------------------------------------------------------
# Mirrors the build_fig1 logic from oneshot_comparison.R, but replaces the
# "Global (whole imaged surface)" denominator with the "Animal cap" region
# (cells with theta <= margin_theta - marg_w, i.e. the tissue on the animal
# side of the margin band).
#
# Produces:
#   * per-frame CSV with densities for {head mesoderm disc, animal cap, global}
#   * per-species binned fold-change CSV (head mesoderm / animal cap)
#   * PDF figure: top view of the animal cap region + density time series +
#     fold change ratio (one panel per normaliser, animal cap vs global).
# =============================================================================

suppressPackageStartupMessages({
  source("renv/activate.R")
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
})

OUT_DIR <- "results/oneshot_comparison"
dir.create(OUT_DIR, showWarnings = FALSE)

MEDAKA_INPUT    <- "data/oriented_medaka_ultrack"
ZEBRAFISH_INPUT <- "data/oriented_zebrafish_ultrack"

MEDAKA_FI  <- 30
ZEB_FI     <- 120
MEDAKA_VOXEL_UM <- 1.05152
ZEB_VOXEL_UM    <- 1.24785

MEDAKA_INGRESSION_FRAME <- 199L
ZEB_INGRESSION_FRAME    <- 40L

# Same animal-cap margin half-widths as oneshot_comparison.R
ZONE_MARG_W_M  <- 5
ZONE_MARG_W_Z  <- 8

# Match the fig1 caps so the timing window is identical to the original
FIG1_ZEB_TCAP_MIN <- 170
FIG1_BASELINE_FROM_INGR_MIN <- -30
FIG1_COMPARE_BIN_MIN <- 2

species_colors <- c("Medaka" = "#E69F00", "Zebrafish" = "#0072B2")
SNAP_NMAX <- 30000

theme_pub <- function(bs = 11) {
  theme_minimal(base_size = bs) +
    theme(
      plot.title    = element_text(face = "bold", size = bs + 1, hjust = 0),
      plot.subtitle = element_text(size = bs - 1, color = "grey40", hjust = 0),
      strip.text    = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey96", color = NA),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

save_pdf <- function(p, name, w = 12, h = 8) {
  path <- file.path(OUT_DIR, name)
  tmp <- tempfile(pattern = "plot_", fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  ggsave(tmp, p, width = w, height = h, device = "pdf")
  if (!file.exists(tmp)) stop("pdf device did not write output: ", tmp)
  ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) stop("failed to overwrite output pdf: ", path)
  cat(sprintf("  saved %s\n", path))
}

banner <- function(x) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", x, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Load tracks
# -----------------------------------------------------------------------------
banner("LOAD INPUT DATA")

sp_m <- fread(file.path(MEDAKA_INPUT, "oriented_tracks_medaka.csv"),
              showProgress = FALSE)
sp_z <- fread(file.path(ZEBRAFISH_INPUT, "oriented_tracks_zebrafish.csv"),
              showProgress = FALSE)

for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z", "RADIAL_DIST",
              "SPHERICAL_DEPTH")) {
  sp_m[[col]] <- sp_m[[col]] * MEDAKA_VOXEL_UM
  sp_z[[col]] <- sp_z[[col]] * ZEB_VOXEL_UM
}

# Same minimum track length filters as the original
ids_m <- sp_m[, .N, by = TRACK_ID][N >= 20, TRACK_ID]
ids_z <- sp_z[, .N, by = TRACK_ID][N >= 5,  TRACK_ID]
sp_m <- sp_m[TRACK_ID %in% ids_m]
sp_z <- sp_z[TRACK_ID %in% ids_z]

sp_m[, `:=`(species = "Medaka",    time_min = FRAME * MEDAKA_FI / 60)]
sp_z[, `:=`(species = "Zebrafish", time_min = FRAME * ZEB_FI    / 60)]

cat(sprintf("  Medaka:    %s spots / %s tracks\n",
            format(nrow(sp_m), big.mark = ","),
            format(uniqueN(sp_m$TRACK_ID), big.mark = ",")))
cat(sprintf("  Zebrafish: %s spots / %s tracks\n",
            format(nrow(sp_z), big.mark = ","),
            format(uniqueN(sp_z$TRACK_ID), big.mark = ",")))

# -----------------------------------------------------------------------------
# Sphere + landmarks + data-driven bulge (copied verbatim from the oneshot)
# -----------------------------------------------------------------------------
read_sphere <- function(dir, vox) {
  d <- fread(file.path(dir, "sphere_params.csv"))
  as.numeric(d[parameter == "radius", value]) * vox
}
R_M <- read_sphere(MEDAKA_INPUT, MEDAKA_VOXEL_UM)
R_Z <- read_sphere(ZEBRAFISH_INPUT, ZEB_VOXEL_UM)

read_landmark <- function(dir, name) {
  d <- fread(file.path(dir, "gastrulation_landmarks.csv"))
  row <- d[landmark == name]
  if (nrow(row) == 0) return(c(theta = NA_real_, phi = NA_real_))
  c(theta = as.numeric(row$theta[1]), phi = as.numeric(row$phi[1]))
}
lm_m_margin <- read_landmark(MEDAKA_INPUT, "margin")
lm_z_margin <- read_landmark(ZEBRAFISH_INPUT, "margin")
MARGIN_M <- lm_m_margin["theta"]
MARGIN_Z <- lm_z_margin["theta"]
lm_m_ingr <- read_landmark(MEDAKA_INPUT, "ingression_center")
lm_z_ingr <- read_landmark(ZEBRAFISH_INPUT, "ingression_center")
cat(sprintf("  Margin theta: medaka %.2f, zebrafish %.2f\n", MARGIN_M, MARGIN_Z))

INGR_DEPTH_PCTILE <- 0.95
detect_bulge <- function(sp, fallback_theta, fallback_phi) {
  t_max <- max(sp$FRAME, na.rm = TRUE)
  late  <- sp[FRAME >= t_max / 2 & is.finite(SPHERICAL_DEPTH)]
  thr <- as.numeric(quantile(late$SPHERICAL_DEPTH, INGR_DEPTH_PCTILE, na.rm = TRUE))
  cand <- late[SPHERICAL_DEPTH >= thr]
  if (nrow(cand) >= 50) {
    ct <- mean(cand$THETA_DEG, na.rm = TRUE)
    cp <- mean(cand$PHI_DEG,   na.rm = TRUE)
    cand[, ang_dist := sqrt((THETA_DEG - ct)^2 + (PHI_DEG - cp)^2)]
    rad <- min(as.numeric(quantile(cand$ang_dist, 0.75, na.rm = TRUE)), 7)
    list(theta = ct, phi = cp, radius = rad, n_cand = nrow(cand),
         depth_thresh = thr, source = "data-driven")
  } else {
    list(theta = fallback_theta, phi = fallback_phi, radius = 10,
         n_cand = nrow(cand), depth_thresh = thr, source = "landmark fallback")
  }
}
BULGE_M <- detect_bulge(sp_m, lm_m_ingr["theta"], lm_m_ingr["phi"])
BULGE_Z <- detect_bulge(sp_z, lm_z_ingr["theta"], lm_z_ingr["phi"])
cat(sprintf("  Bulge medaka:    theta=%.1f phi=%.1f r=%.1f (n=%d)\n",
            BULGE_M$theta, BULGE_M$phi, BULGE_M$radius, BULGE_M$n_cand))
cat(sprintf("  Bulge zebrafish: theta=%.1f phi=%.1f r=%.1f (n=%d)\n",
            BULGE_Z$theta, BULGE_Z$phi, BULGE_Z$radius, BULGE_Z$n_cand))

# -----------------------------------------------------------------------------
# Animal cap definition: theta <= margin_theta - marg_w
# Spherical cap area (from theta=0 to theta=theta_cap):  2*pi*R^2*(1 - cos(t))
# -----------------------------------------------------------------------------
THETA_CAP_M <- MARGIN_M - ZONE_MARG_W_M
THETA_CAP_Z <- MARGIN_Z - ZONE_MARG_W_Z

animal_cap_area <- function(R, theta_cap_deg) {
  t_rad <- theta_cap_deg * pi / 180
  2 * pi * R^2 * (1 - cos(t_rad))
}
A_ANIM_M <- animal_cap_area(R_M, THETA_CAP_M)
A_ANIM_Z <- animal_cap_area(R_Z, THETA_CAP_Z)
cat(sprintf("  Animal cap boundary: medaka theta<=%.1f (A=%.1f x10^3 um^2), zebrafish theta<=%.1f (A=%.1f x10^3 um^2)\n",
            THETA_CAP_M, A_ANIM_M/1000, THETA_CAP_Z, A_ANIM_Z/1000))

# Head-mesoderm disc area (same as oneshot: ellipse with semi-axes R*r_rad in
# theta and R*sin(theta0)*r_rad in phi)
ellipse_disc_area <- function(R, b)
  pi * R^2 * (b$radius * pi / 180)^2 * sin(b$theta * pi / 180)
A_BULGE_M  <- ellipse_disc_area(R_M, BULGE_M)
A_BULGE_Z  <- ellipse_disc_area(R_Z, BULGE_Z)

# Global = imaged surface area (used for the side-by-side panel)
imaged_area_um2 <- function(dt, R) {
  th_q <- as.numeric(quantile(dt$THETA_DEG, c(0.02, 0.98), na.rm = TRUE)) * pi / 180
  ph_q <- as.numeric(quantile(dt$PHI_DEG,   c(0.02, 0.98), na.rm = TRUE)) * pi / 180
  R^2 * (th_q[2] - th_q[1]) * sin(mean(th_q)) * (ph_q[2] - ph_q[1])
}
A_GLOBAL_M <- imaged_area_um2(sp_m, R_M)
A_GLOBAL_Z <- imaged_area_um2(sp_z, R_Z)

# -----------------------------------------------------------------------------
# Selectors + per-frame density (depth-integrated, like the original)
# -----------------------------------------------------------------------------
in_disc <- function(sp, b) sp[!is.na(PHI_DEG) &
  sqrt((THETA_DEG - b$theta)^2 + (PHI_DEG - b$phi)^2) < b$radius]

in_animal_cap <- function(sp, theta_cap) {
  sp[!is.na(THETA_DEG) & THETA_DEG <= theta_cap]
}

per_frame_density <- function(cells, A, fi, sp_lbl, reg_lbl) {
  cells[, .(n = .N), by = FRAME][,
        .(species = sp_lbl, region = reg_lbl, FRAME,
          time_min = FRAME * fi / 60,
          n_cells = n,
          area_um2 = A,
          density_per_1000um2 = n / A * 1000)]
}

FIG1_TIME_META <- data.table(
  species = c("Medaka", "Zebrafish"),
  ingression_frame = c(MEDAKA_INGRESSION_FRAME, ZEB_INGRESSION_FRAME),
  fi_sec = c(MEDAKA_FI, ZEB_FI))
FIG1_TIME_META[, ingression_time_min := ingression_frame * fi_sec / 60]

REGION_LVL <- c("Head mesoderm disc", "Animal cap", "Global (whole imaged surface)")
REGION_COL <- c("Head mesoderm disc"               = "#E31A1C",
                "Animal cap"                       = "#3B7DD8",
                "Global (whole imaged surface)"   = "grey40")

# -----------------------------------------------------------------------------
# Build per-frame densities for both species, both versions
# -----------------------------------------------------------------------------
banner("BUILD DENSITY TIME SERIES")

# Full version (medaka = full, zebrafish = t<=170 min)
sp_m_t <- copy(sp_m)
sp_z_t <- sp_z[time_min <= FIG1_ZEB_TCAP_MIN]

d_ts <- rbind(
  per_frame_density(in_disc(sp_m_t, BULGE_M),        A_BULGE_M,  MEDAKA_FI, "Medaka",    "Head mesoderm disc"),
  per_frame_density(in_disc(sp_z_t, BULGE_Z),        A_BULGE_Z,  ZEB_FI,    "Zebrafish", "Head mesoderm disc"),
  per_frame_density(in_animal_cap(sp_m_t, THETA_CAP_M), A_ANIM_M, MEDAKA_FI, "Medaka",    "Animal cap"),
  per_frame_density(in_animal_cap(sp_z_t, THETA_CAP_Z), A_ANIM_Z, ZEB_FI,    "Zebrafish", "Animal cap"),
  per_frame_density(sp_m_t,                           A_GLOBAL_M, MEDAKA_FI, "Medaka",    "Global (whole imaged surface)"),
  per_frame_density(sp_z_t,                           A_GLOBAL_Z, ZEB_FI,    "Zebrafish", "Global (whole imaged surface)"))

d_ts <- merge(d_ts,
              FIG1_TIME_META[, .(species, ingression_frame, ingression_time_min)],
              by = "species", all.x = TRUE)
d_ts[, time_from_ingression_min := time_min - ingression_time_min]
d_ts[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
d_ts[, region := factor(region, levels = REGION_LVL)]

# QC: drop frames where global density < 50% of species median
global_qc <- d_ts[region == "Global (whole imaged surface)",
                  .(thr = 0.5 * median(density_per_1000um2)), by = species]
bad <- d_ts[region == "Global (whole imaged surface)"][
  global_qc, on = "species"][density_per_1000um2 < thr, .(species, FRAME)]
if (nrow(bad)) {
  cat(sprintf("  dropped %d broken frames (global < 50%% of median)\n", nrow(bad)))
  d_ts <- d_ts[!bad, on = c("species", "FRAME")]
}

fwrite(d_ts, file.path(OUT_DIR, "head_meso_vs_animal_pole_per_frame.csv"))

# -----------------------------------------------------------------------------
# Fold change vs each normaliser:  (head mesoderm density / head mesoderm base)
#                                    / (NORM density / NORM base)
# Baseline window: [-30, 0) min from ingression (same as the original)
# -----------------------------------------------------------------------------
baseline_dt <- d_ts[, {
  dens <- density_per_1000um2[
    time_from_ingression_min >= FIG1_BASELINE_FROM_INGR_MIN &
    time_from_ingression_min < 0
  ]
  if (!length(dens)) dens <- density_per_1000um2[time_from_ingression_min < 0]
  if (!length(dens)) dens <- head(density_per_1000um2, 5)
  .(baseline_density = mean(dens, na.rm = TRUE))
}, by = .(species, region)]

fold_ratio <- function(d_ts_in, baseline_dt_in, numerator_region, denom_region) {
  num <- d_ts_in[region == numerator_region,
                 .(species, FRAME, time_min, time_from_ingression_min,
                   num_density = density_per_1000um2)]
  num_b <- baseline_dt_in[region == numerator_region,
                          .(species, num_base = baseline_density)]
  den <- d_ts_in[region == denom_region,
                 .(species, FRAME, den_density = density_per_1000um2)]
  den_b <- baseline_dt_in[region == denom_region,
                          .(species, den_base = baseline_density)]
  m <- merge(num, num_b, by = "species")
  m <- merge(m, den, by = c("species", "FRAME"))
  m <- merge(m, den_b, by = "species")
  m[, num_fold := num_density / num_base]
  m[, den_fold := den_density / den_base]
  m[, fold_ratio := num_fold / den_fold]
  m[, `:=`(numerator_region = numerator_region,
           denom_region = denom_region)]
  m[]
}

fr_anim <- fold_ratio(d_ts, baseline_dt,
                      numerator_region = "Head mesoderm disc",
                      denom_region     = "Animal cap")
fr_glob <- fold_ratio(d_ts, baseline_dt,
                      numerator_region = "Head mesoderm disc",
                      denom_region     = "Global (whole imaged surface)")

# Combined (long) for plotting
fr_long <- rbindlist(list(
  fr_anim[, .(species, FRAME, time_min, time_from_ingression_min,
              fold_ratio, numerator_region, denom_region, normaliser = "Animal cap")],
  fr_glob[, .(species, FRAME, time_min, time_from_ingression_min,
              fold_ratio, numerator_region, denom_region,
              normaliser = "Global (whole imaged surface)")]),
  use.names = TRUE)
fr_long[, normaliser := factor(normaliser,
                               levels = c("Animal cap",
                                          "Global (whole imaged surface)"))]

# Binned at 2-min from-ingression steps (matches original)
fr_binned <- fr_long[, .(
  fold_ratio = mean(fold_ratio, na.rm = TRUE),
  n_obs = .N),
  by = .(species, normaliser,
         time_plot_min = floor(time_from_ingression_min / FIG1_COMPARE_BIN_MIN) *
           FIG1_COMPARE_BIN_MIN)]

ratio_summary <- fr_binned[, .(
  baseline_window = sprintf("[%d, 0) min", FIG1_BASELINE_FROM_INGR_MIN),
  peak_ratio = round(max(fold_ratio, na.rm = TRUE), 3),
  peak_time_from_ingression = time_plot_min[which.max(fold_ratio)]
), by = .(species, normaliser)]

cat("Fold-change vs each normaliser (peak, baseline window = [-30, 0) min):\n")
print(ratio_summary)
fwrite(fr_long,   file.path(OUT_DIR, "head_meso_vs_animal_pole_fold_long.csv"))
fwrite(fr_binned, file.path(OUT_DIR, "head_meso_vs_animal_pole_fold_binned.csv"))
fwrite(ratio_summary, file.path(OUT_DIR, "head_meso_vs_animal_pole_summary.csv"))

# -----------------------------------------------------------------------------
# Per-region fold change (vs each region's own baseline) -- no division.
# This shows the underlying trajectories that drive the ratios in C / D.
# -----------------------------------------------------------------------------
fc_long <- d_ts[, .(species, FRAME, time_min, time_from_ingression_min,
                    region, density_per_1000um2)][
  baseline_dt, on = .(species, region)][
    , fold := density_per_1000um2 / baseline_density][
      , .(species, FRAME, time_min, time_from_ingression_min,
          region, fold)]
fc_binned <- fc_long[, .(fold = mean(fold, na.rm = TRUE), n_obs = .N),
                     by = .(species, region,
                            time_plot_min = floor(time_from_ingression_min /
                                                    FIG1_COMPARE_BIN_MIN) *
                              FIG1_COMPARE_BIN_MIN)]
fwrite(fc_long,   file.path(OUT_DIR, "head_meso_vs_animal_pole_fc_long.csv"))
fwrite(fc_binned, file.path(OUT_DIR, "head_meso_vs_animal_pole_fc_binned.csv"))

# -----------------------------------------------------------------------------
# Absolute density (nuclei / 1000 um^2) of head mesoderm disc and animal cap,
# both species, on the same axes.  No baseline, no division -- direct
# comparison in the same physical units.
# -----------------------------------------------------------------------------
abs_long <- d_ts[region %in% c("Head mesoderm disc", "Animal cap"),
                 .(species, region, FRAME, time_min, time_from_ingression_min,
                   density_per_1000um2)]
abs_binned <- abs_long[, .(density_per_1000um2 = mean(density_per_1000um2,
                                                      na.rm = TRUE),
                            n_obs = .N),
                        by = .(species, region,
                               time_plot_min = floor(time_from_ingression_min /
                                                       FIG1_COMPARE_BIN_MIN) *
                                 FIG1_COMPARE_BIN_MIN)]
fwrite(abs_long,   file.path(OUT_DIR,
                              "head_meso_vs_animal_pole_abs_long.csv"))
fwrite(abs_binned, file.path(OUT_DIR,
                              "head_meso_vs_animal_pole_abs_binned.csv"))

# -----------------------------------------------------------------------------
# FIGURE
#   A. Top view of embryo surface with the animal-cap region shaded
#      (per species) plus the head-mesoderm disc outlined.
#   B. Density time series for {head mesoderm, animal cap, global}.
#   C. Fold change: head mesoderm / animal cap  (the new normaliser)
#   D. Fold change: head mesoderm / global       (for direct comparison)
#   E. Per-region fold change over baseline (no division) — underlying
#      trajectories of all three regions side by side.
# -----------------------------------------------------------------------------
banner("BUILD FIGURE")

# Panel A: snapshot + animal cap + bulge disc
late_snap <- function(sp, fi, t_cap, w = 10) {
  d <- sp[time_min <= t_cap]
  tmax <- max(d$time_min, na.rm = TRUE)
  d[time_min >= tmax - w]
}
sample_n <- function(d, n) if (nrow(d) > n) d[sample(.N, n)] else d

disc_outline <- function(b, sp_lbl, n = 80) {
  ang <- seq(0, 2 * pi, length.out = n)
  data.table(species = sp_lbl,
             phi   = b$phi   + b$radius * cos(ang),
             theta = b$theta + b$radius * sin(ang))
}
anim_band_rect <- data.table(
  species = factor(c("Medaka", "Zebrafish"), levels = c("Medaka", "Zebrafish")),
  ymin    = c(-Inf, -Inf),
  ymax    = c(THETA_CAP_M, THETA_CAP_Z))

snap_m <- late_snap(sp_m, MEDAKA_FI, Inf)
snap_z <- late_snap(sp_z, ZEB_FI,    FIG1_ZEB_TCAP_MIN)
tA_m <- max(snap_m$time_min); tA_z <- max(snap_z$time_min)
snap <- rbind(
  sample_n(snap_m, SNAP_NMAX)[, .(species = "Medaka",    PHI_DEG, THETA_DEG, SPHERICAL_DEPTH)],
  sample_n(snap_z, SNAP_NMAX)[, .(species = "Zebrafish", PHI_DEG, THETA_DEG, SPHERICAL_DEPTH)])
snap[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

outlines <- rbind(disc_outline(BULGE_M, "Medaka"),
                  disc_outline(BULGE_Z, "Zebrafish"))
outlines[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

pA <- ggplot(snap, aes(PHI_DEG, THETA_DEG)) +
  geom_rect(data = anim_band_rect, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = REGION_COL[["Animal cap"]], alpha = 0.18) +
  geom_point(aes(color = SPHERICAL_DEPTH), size = 0.25, alpha = 0.55) +
  geom_polygon(data = outlines, inherit.aes = FALSE,
               aes(phi, theta), fill = NA, color = "#E31A1C", linewidth = 0.9) +
  facet_wrap(~ species, scales = "free") +
  scale_y_reverse() +
  scale_color_viridis_c(option = "inferno", name = expression("depth ("*mu*"m)")) +
  labs(title = sprintf("A. Top view of the embryo surface (snapshot: medaka %.0f-%.0f min, zebrafish %.0f-%.0f min)",
                       tA_m - 10, tA_m, tA_z - 10, tA_z),
       subtitle = sprintf("Blue band = animal cap (theta <= margin - %d deg medaka / %d deg zebrafish). Red = head-mesoderm / bulge disc.",
                          ZONE_MARG_W_M, ZONE_MARG_W_Z),
       x = expression(varphi*" (deg)"),
       y = expression(theta*" (deg, animal pole up)")) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"))

# Panel B: density vs real time
pB <- ggplot(d_ts, aes(time_min, density_per_1000um2, color = region)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, linewidth = 1) +
  facet_wrap(~ species, nrow = 1, scales = "free_x") +
  scale_color_manual(values = REGION_COL, name = NULL) +
  labs(title = "B. Regional nuclei density trajectories over real time",
       subtitle = "Head-mesoderm disc (red), animal cap (blue), global imaged surface (grey). Same depth-integrated definition as Fig 1.",
       x = "time (min)",
       y = expression("nuclei / 1000 "*mu*"m"^2)) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        legend.position = "bottom")

# Panels C / D: aligned fold change
plot_panel_cd <- function(df, normaliser_lbl, title_lbl, subtitle_lbl) {
  d <- df[normaliser == normaliser_lbl]
  ggplot(d, aes(time_plot_min, fold_ratio, color = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey65") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey55") +
    geom_line(linewidth = 0.9, alpha = 0.95) +
    geom_point(size = 1.3, alpha = 0.85) +
    scale_color_manual(values = species_colors, name = NULL) +
    labs(title = title_lbl, subtitle = subtitle_lbl,
         x = "time from ingression (min)",
         y = "head-mesoderm fold change / normaliser fold change") +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"))
}

pC <- plot_panel_cd(
  fr_binned, "Animal cap", "C. Head mesoderm disc / Animal cap",
  sprintf("Local fold change divided by animal-cap fold change.  Animal cap = theta <= margin - %d deg (medaka) / %d deg (zebrafish).  2-min aligned bins.  Baseline = [%d, 0) min.",
          ZONE_MARG_W_M, ZONE_MARG_W_Z, FIG1_BASELINE_FROM_INGR_MIN))

pD <- plot_panel_cd(
  fr_binned, "Global (whole imaged surface)",
  "D. Head mesoderm disc / Global (reproduces Fig 1D for direct comparison)",
  sprintf("Same metric using the global imaged surface as normaliser.  Baseline = [%d, 0) min.",
          FIG1_BASELINE_FROM_INGR_MIN))

# Panel E: ABSOLUTE density (no baseline, no division) -- head mesoderm disc,
# both species, on the same axes, coloured by species.
pE <- ggplot(abs_binned[region == "Head mesoderm disc"],
             aes(time_plot_min, density_per_1000um2, color = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  scale_color_manual(values = species_colors, name = NULL) +
  labs(title = "E. Absolute head-mesoderm disc density (no baseline, no division)",
       subtitle = "Raw nuclei per 1000 um^2, depth-integrated over the head-mesoderm disc.  Both species on the same y-axis.  2-min aligned bins from ingression.",
       x = "time from ingression (min)",
       y = expression("nuclei / 1000 "*mu*"m"^2),
       color = NULL) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 9, color = "grey25"),
        legend.position = "bottom")

fig <- (pA / pB / (pC | pD) / pE) +
  plot_layout(heights = c(1.2, 1.1, 1.1, 1.1)) +
  plot_annotation(
    title = "Head-mesoderm / bulge density normalised by the ANIMAL POLE (animal cap)",
    subtitle = sprintf("Each curve is the local density fold change from a [-%d, 0) min baseline, divided by the normaliser's own fold change.  >1 = disc-specific enrichment beyond the normaliser's own densification.  0 min = ingression (medaka frame %d = %.1f min; zebrafish frame %d = %.1f min).",
                       -FIG1_BASELINE_FROM_INGR_MIN,
                       MEDAKA_INGRESSION_FRAME,
                       FIG1_TIME_META[species == "Medaka", ingression_time_min],
                       ZEB_INGRESSION_FRAME,
                       FIG1_TIME_META[species == "Zebrafish", ingression_time_min]),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, color = "grey25")))

save_pdf(fig, "01c_head_meso_vs_animal_pole.pdf", w = 16, h = 24)

# Peak summary
cat("\nPeak fold change vs each normaliser:\n")
print(ratio_summary)
banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
