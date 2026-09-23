suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

ZEB_DIR  <- "data/oriented_zebrafish_ultrack"
ZEB_TRACK_FILE <- "oriented_tracks_zebrafish.csv"
ZEB_VOX  <- 1.24785
ZEB_FI   <- 120  # s

sp <- fread(file.path(ZEB_DIR, ZEB_TRACK_FILE))
# unify column names if needed
setnames(sp, old = grep("^THETA", names(sp), value = TRUE)[1], new = "THETA_DEG", skip_absent = TRUE)
setnames(sp, old = grep("^PHI",   names(sp), value = TRUE)[1], new = "PHI_DEG",   skip_absent = TRUE)
setnames(sp, old = grep("SPHERICAL_DEPTH", names(sp), value = TRUE)[1], new = "SPHERICAL_DEPTH", skip_absent = TRUE)
if (!"FRAME" %in% names(sp)) setnames(sp, grep("FRAME", names(sp), ignore.case=TRUE, value=TRUE)[1], "FRAME")

sp[, time_min := FRAME * ZEB_FI / 60]
sp_dir <- ZEB_DIR
R <- as.numeric(fread(file.path(sp_dir, "sphere_params.csv"))[parameter=="radius", value]) * ZEB_VOX
lm  <- fread(file.path(sp_dir, "gastrulation_landmarks.csv"))
M_TH <- lm[landmark=="margin", theta][1]
cat(sprintf("Zebrafish: sphere R=%.1f um, MARGIN_Z=%.2f deg\n", R, M_TH))
cat(sprintf("Spots: %s, frames %d-%d, %.0f min total\n",
            format(nrow(sp), big.mark=","), min(sp$FRAME), max(sp$FRAME), max(sp$time_min)))

# theta range / occupancy
cat("\nTheta percentiles across all spots:\n")
print(round(quantile(sp$THETA_DEG, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm=TRUE), 1))

# At t=0 vs t=170 vs end: where do cells live in theta? what is depth profile?
times <- c(0, 60, 120, 170, max(sp$time_min, na.rm=TRUE))
cat("\nTheta distribution at selected times (+-5 min window):\n")
for (t in times) {
  d <- sp[time_min >= t-5 & time_min <= t+5]
  q <- round(quantile(d$THETA_DEG, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), 1)
  cat(sprintf("  t=%6.1f min  n=%6d  theta p05/p25/p50/p75/p95 = %s\n",
              t, nrow(d), paste(q, collapse=" / ")))
}

# Margin slab definition
W <- 8
cat(sprintf("\nMargin slab: theta in [%.1f, %.1f] (half-width %d deg)\n", M_TH-W, M_TH+W, W))
# How many cells lie INSIDE that slab at each time?
sp[, in_margin := THETA_DEG >= (M_TH-W) & THETA_DEG <= (M_TH+W)]
ct <- sp[, .(n = .N, n_margin = sum(in_margin),
             frac_margin = mean(in_margin),
             depth_p10 = quantile(SPHERICAL_DEPTH, 0.10, na.rm=TRUE),
             depth_p50 = quantile(SPHERICAL_DEPTH, 0.50, na.rm=TRUE),
             depth_p90 = quantile(SPHERICAL_DEPTH, 0.90, na.rm=TRUE),
             depth_range = quantile(SPHERICAL_DEPTH, 0.90, na.rm=TRUE)-quantile(SPHERICAL_DEPTH, 0.10, na.rm=TRUE)),
         by = .(time_bin = floor(time_min/20)*20)][order(time_bin)]
cat("\nGlobal occupancy and margin-slab statistics by 20-min bin:\n")
print(ct, nrows=50)

# Within-margin only: depth percentiles per time bin
mg <- sp[in_margin == TRUE,
         .(n = .N,
           depth_min = min(SPHERICAL_DEPTH, na.rm=TRUE),
           depth_p10 = quantile(SPHERICAL_DEPTH, 0.10, na.rm=TRUE),
           depth_p25 = quantile(SPHERICAL_DEPTH, 0.25, na.rm=TRUE),
           depth_p50 = quantile(SPHERICAL_DEPTH, 0.50, na.rm=TRUE),
           depth_p75 = quantile(SPHERICAL_DEPTH, 0.75, na.rm=TRUE),
           depth_p90 = quantile(SPHERICAL_DEPTH, 0.90, na.rm=TRUE),
           depth_max = max(SPHERICAL_DEPTH, na.rm=TRUE),
           range_p10_p90 = quantile(SPHERICAL_DEPTH,0.90,na.rm=TRUE)-quantile(SPHERICAL_DEPTH,0.10,na.rm=TRUE)),
         by = .(time_bin = floor(time_min/10)*10)][order(time_bin)]
cat("\nMargin-slab depth distribution by 10-min bin:\n")
print(mg, nrows=60)

# Find which theta band is actually the THICKEST band in the embryo at each time --
# this tells whether the chosen MARGIN_Z is correctly placed.
BAND <- 4  # half width
sp[, theta_band := round(THETA_DEG/BAND)*BAND]
band_thick <- sp[is.finite(SPHERICAL_DEPTH),
   .(n = .N,
     range_p10_p90 = quantile(SPHERICAL_DEPTH,0.90,na.rm=TRUE)-quantile(SPHERICAL_DEPTH,0.10,na.rm=TRUE),
     p90 = quantile(SPHERICAL_DEPTH, 0.90, na.rm=TRUE)),
   by = .(theta_band, time_bin = floor(time_min/30)*30)][n >= 30]

cat(sprintf("\nThickest theta-band (range p10-p90) per 30-min window  -- MARGIN_Z=%.1f deg:\n", M_TH))
for (tt in sort(unique(band_thick$time_bin))) {
  d <- band_thick[time_bin == tt][order(-range_p10_p90)][1:5]
  cat(sprintf("  t=%3d min  top-5 thick bands (theta -> range p10-p90 um):\n", tt))
  for (i in seq_len(nrow(d))) {
    cat(sprintf("    theta %5.1f  range %5.1f um  (n=%d)\n",
                d$theta_band[i], d$range_p10_p90[i], d$n[i]))
  }
}

# Plot: heat-map theta vs time, coloured by range_p10_p90
band_thick[, theta_band := as.numeric(theta_band)]
pp <- ggplot(band_thick, aes(time_bin, theta_band, fill = range_p10_p90)) +
  geom_tile() +
  geom_hline(yintercept = M_TH, color = "white", linewidth = 0.6) +
  geom_hline(yintercept = c(M_TH-W, M_TH+W), color = "white", linetype="dashed", linewidth = 0.4) +
  scale_y_reverse() +
  scale_fill_viridis_c(option="inferno", name = "p90-p10 depth (um)") +
  labs(title = sprintf("Zebrafish — depth range (p90-p10) per theta band (%d deg) over time", BAND),
       subtitle = sprintf("White solid line = MARGIN_Z = %.1f deg.  Dashed = margin slab +/-%d deg.", M_TH, W),
       x = "time (min)", y = "theta (deg, animal pole at 0)") +
  theme_minimal()

ggsave("results/oneshot_comparison/zeb_margin_check.pdf", pp, width=10, height=6)
cat("\nSaved results/oneshot_comparison/zeb_margin_check.pdf\n")
