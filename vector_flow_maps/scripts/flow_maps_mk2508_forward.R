# =============================================================================
# FLOW MAPS - mk2508  (INTERNALISATION anchors only, 30 min FORWARD window)
# =============================================================================
# Adapted from flow_maps_mk2508_stephane.R:
#   * Uses ONLY the internalisation reference frames (3 landmarks; the
#     contraction-anchored panel set is dropped entirely).
#   * Each panel averages the WIN_MIN minutes STARTING +1 frame AFTER each
#     landmark (forward window = [ref_frame + 1, ref_frame + WIN_FRAMES]).
#     This is the temporal mirror of the original "30 min before the
#     reference" panel direction.
#
# Reference frames (absolute frame units, supplied by user; FI = 30 s/frame
# = 0.5 min/frame):
#   internalisation_minus_1h = 40L     (1 h before internalisation)
#   internalisation_start    = 160L    (start of internalisation)
#   internalisation_plus_1h  = 280L    (1 h after internalisation)
#
# Forward 30 min (= 60 frames) windows for WIN_MIN = 30:
#   ref=40   -> window frames 41-100
#   ref=160  -> window frames 161-220
#   ref=280  -> window frames 281-340
#
# Inputs:
#   inputs/tracks_for_flow_maps/tracks_mk_2508_stephane_oriented_filtered.csv
#
# Outputs (in outputs/mk2508/forward_internalisation/):
#   04a_flow_<WIN_MIN>min_after_internalization_mk2508.pdf   (3 panels)
#   04c_majority_movement_mk2508_forward.pdf                  (majority strip)
#   flow_fields.csv
#   panel_metadata.csv
#   panel_majority_movement.csv
#   summary_references.txt
#
# Run from the vector_flow_maps/ root:
#   Rscript scripts/flow_maps_mk2508_forward.R 30
# =============================================================================

suppressPackageStartupMessages({
  source("renv/activate.R")
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

OUT_DIR  <- file.path("outputs", "mk2508", "forward_internalisation")
IN_CSV   <- file.path("inputs", "tracks_for_flow_maps",
                      "tracks_mk_2508_stephane_oriented_filtered.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

VOXEL_UM <- 0.64
FI_SEC   <- 30
FI_MIN   <- FI_SEC / 60

# Forward-window length (in minutes), can be overridden via CLI:
#   Rscript flow_maps_mk2508_forward.R 30
WIN_MIN    <- as.integer(commandArgs(trailingOnly = TRUE)[1] %||% 30)
WIN_FRAMES <- as.integer(WIN_MIN / FI_MIN)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# Internalisation anchors ONLY (no contraction anchors).
INTERNALIZATION_REFERENCE_FRAMES <- list(
  internalization_minus_1h = 40L,
  internalization_start    = 160L,
  internalization_plus_1h  = 280L
)
INTERNALIZATION_LABELS <- c(
  "1 h before internalization",
  "Start of internalization",
  "1 h after internalization"
)

# Flow-field tuning (identical to original)
FLOW_BIN_UM       <- 30
FLOW_MIN_N        <- 10
ARROW_SCALE       <- 60
ARROW_HEAD        <- 0.14
ARROW_LW          <- 0.7
VORT_SWIRL_THRESH <- 0.3
SMOOTH_K          <- 5L
VEL_LAG           <- 4L     # 2-minute step
MIN_FRAMES        <- 20L

flow_cols <- c("Upward (AP)" = "#2166AC", "Circular" = "#1A9850",
               "Downward (VP)" = "#B2182B")

theme_pub <- function(bs = 10) {
  theme_minimal(base_size = bs) +
    theme(
      plot.title    = element_text(face = "bold", size = bs + 1, hjust = 0),
      plot.subtitle = element_text(size = bs - 1, color = "grey40", hjust = 0),
      strip.text    = element_text(face = "bold", size = 8),
      strip.background = element_rect(fill = "grey96", color = NA),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

save_pdf <- function(p, name, w = 18, h) {
  path <- file.path(OUT_DIR, name)
  tmp  <- tempfile(pattern = "plot_", fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  ggsave(tmp, p, width = w, height = h, device = "pdf")
  if (!file.exists(tmp)) stop("pdf device did not write output: ", tmp)
  ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) stop("failed to overwrite output pdf: ", path)
  cat(sprintf("  saved %s\n", path))
}

ascii_safe <- function(s) {
  s <- gsub("\u0394", "d", s, fixed = TRUE)
  s <- gsub("\u2212", "-", s, fixed = TRUE)
  s <- gsub("\u2014", "-", s, fixed = TRUE)
  s <- gsub("\u2013", "-", s, fixed = TRUE)
  s <- gsub("\u00B1", "+/-", s, fixed = TRUE)
  iconv(s, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "?")
}

banner <- function(x) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", x, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Load + orient
# -----------------------------------------------------------------------------

banner("LOAD ORIENTED + FILTERED TRACKS")

sp <- fread(IN_CSV, showProgress = FALSE)
cat(sprintf("  loaded %s rows, %d columns\n", format(nrow(sp), big.mark = ","),
            ncol(sp)))
cat(sprintf("  columns: %s\n", paste(names(sp), collapse = ", ")))

if ("FILTERED_OUT" %in% names(sp)) {
  drop_mask <- tolower(as.character(sp$FILTERED_OUT)) == "true"
  n_in   <- nrow(sp)
  n_drop <- sum(drop_mask)
  n_keep <- n_in - n_drop
  cat(sprintf("  honouring FILTERED_OUT flag: %s in, %s kept, %s dropped\n",
              format(n_in, big.mark = ","),
              format(n_keep, big.mark = ","),
              format(n_drop, big.mark = ",")))
  sp <- sp[!drop_mask]
  sp[, FILTERED_OUT := NULL]
}

for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z",
              "RADIAL_DIST", "SPHERICAL_DEPTH")) {
  sp[[col]] <- sp[[col]] * VOXEL_UM
}

track_n <- sp[, .N, by = TRACK_ID]
keep_ids <- track_n[N >= MIN_FRAMES, TRACK_ID]
sp <- sp[TRACK_ID %in% keep_ids]

sp[, time_min := FRAME * FI_MIN]
cat(sprintf("  after short-track filter: %s spots / %s tracks (frames %d-%d)\n",
            format(nrow(sp), big.mark = ","),
            format(uniqueN(sp$TRACK_ID), big.mark = ","),
            min(sp$FRAME), max(sp$FRAME)))

# -----------------------------------------------------------------------------
# Velocity / per-step metrics
# -----------------------------------------------------------------------------

banner("COMPUTE PER-STEP VELOCITIES")

compute_vel <- function(df, fi_min, vel_lag, smooth_k) {
  dt <- copy(df); setkey(dt, TRACK_ID, FRAME)
  for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z",
                "RADIAL_DIST", "THETA_DEG", "PHI_DEG", "SPHERICAL_DEPTH")) {
    sm <- paste0(col, "_SM")
    dt[, (sm) := frollmean(get(col), n = smooth_k, align = "center"),
       by = TRACK_ID]
    dt[is.na(get(sm)), (sm) := get(col)]
  }
  dt[, `:=`(
    dx   = POSITION_X_SM - shift(POSITION_X_SM, vel_lag),
    dy   = POSITION_Y_SM - shift(POSITION_Y_SM, vel_lag),
    dz   = POSITION_Z_SM - shift(POSITION_Z_SM, vel_lag),
    dt_f = FRAME - shift(FRAME, vel_lag)
  ), by = TRACK_ID]
  dt[, disp_3d := sqrt(dx^2 + dy^2 + dz^2)]
  dt[, `:=`(
    inst_speed  = disp_3d / (dt_f * fi_min),
    vx_um_min   = dx / (dt_f * fi_min),
    vy_um_min   = dy / (dt_f * fi_min)
  )]
  dt[!is.na(dt_f) & dt_f == vel_lag]
}

vel <- compute_vel(sp, FI_MIN, VEL_LAG, SMOOTH_K)
cat(sprintf("  per-step rows: %s (lag = %d frames = %d min)  frames %d-%d\n",
            format(nrow(vel), big.mark = ","), VEL_LAG, VEL_LAG * FI_MIN,
            min(vel$FRAME), max(vel$FRAME)))

# -----------------------------------------------------------------------------
# FORWARD-window flow-field construction
# -----------------------------------------------------------------------------

banner("BUILD FORWARD FLOW FIELDS (30 min AFTER each reference)")

build_window_flow <- function(vel, ref_frame, win_frames = WIN_FRAMES,
                              bin_um = FLOW_BIN_UM, min_n = FLOW_MIN_N,
                              swirl_thr = VORT_SWIRL_THRESH) {
  # FORWARD window: the WIN_MIN minutes AFTER the reference frame, not
  # including the reference itself.  frame_lo = ref_frame + 1,
  # frame_hi = ref_frame + win_frames (clipped to the available recording).
  frame_lo     <- ref_frame + 1L
  frame_hi_req <- ref_frame + win_frames
  frame_hi     <- min(frame_hi_req, max(vel$FRAME))

  sub <- vel[FRAME >= frame_lo & FRAME <= frame_hi]
  if (nrow(sub) < 200) return(NULL)

  bs <- bin_um
  fl <- sub[, .(
    mean_vx    = mean(vx_um_min, na.rm = TRUE),
    mean_vy    = mean(vy_um_min, na.rm = TRUE),
    mean_speed = mean(inst_speed, na.rm = TRUE),
    n          = .N
  ), by = .(
    x_bin = floor(POSITION_X / bs) * bs + bs / 2,
    y_bin = floor(POSITION_Y / bs) * bs + bs / 2
  )][n >= min_n]
  if (nrow(fl) < 5) return(NULL)

  setkey(fl, x_bin, y_bin)
  svx <- numeric(nrow(fl)); svy <- numeric(nrow(fl))
  for (i in seq_len(nrow(fl))) {
    xb <- fl$x_bin[i]; yb <- fl$y_bin[i]
    vxs <- c(); vys <- c()
    for (dx_ in c(-bs, 0, bs)) for (dy_ in c(-bs, 0, bs)) {
      nb <- fl[.(xb + dx_, yb + dy_)]
      if (nrow(nb) == 1L) { vxs <- c(vxs, nb$mean_vx); vys <- c(vys, nb$mean_vy) }
    }
    svx[i] <- mean(vxs); svy[i] <- mean(vys)
  }
  fl[, c("svx", "svy") := .(svx, svy)]

  vort <- rep(NA_real_, nrow(fl))
  for (i in seq_len(nrow(fl))) {
    xb <- fl$x_bin[i]; yb <- fl$y_bin[i]
    e <- fl[.(xb + bs, yb)]; w <- fl[.(xb - bs, yb)]
    n_nb <- fl[.(xb, yb - bs)]; s_nb <- fl[.(xb, yb + bs)]
    dvy_dx <- if (nrow(e) == 1 && nrow(w) == 1)
                (e$svy - w$svy) / (2 * bs) else NA_real_
    dvx_dy <- if (nrow(n_nb) == 1 && nrow(s_nb) == 1)
                (s_nb$svx - n_nb$svx) / (2 * bs) else NA_real_
    if (!is.na(dvy_dx) && !is.na(dvx_dy)) vort[i] <- dvy_dx - dvx_dy
  }
  fl[, vorticity := vort]
  fl[, speed_2d := sqrt(svx^2 + svy^2)]
  fl[, swirl := abs(vorticity) * bs / pmax(speed_2d, 1e-8)]
  fl[, move_type := fifelse(
    !is.na(swirl) & swirl > swirl_thr, "Circular",
    fifelse(mean_vy > 0, "Downward (VP)", "Upward (AP)"))]
  fl[, move_type := factor(move_type,
      levels = c("Upward (AP)", "Circular", "Downward (VP)"))]
  fl[, frame_lo       := frame_lo]
  fl[, frame_hi       := frame_hi]
  fl[, ref_frame      := ref_frame]
  fl[, window_clipped := frame_hi < frame_hi_req]
  fl[]
}

# Build one flow field per internalisation anchor.
flow_panels <- list()
panel_meta  <- list()

INTERN_FRAMES <- as.integer(unlist(INTERNALIZATION_REFERENCE_FRAMES))
MAX_FRAME     <- max(vel$FRAME)

for (i in seq_along(INTERN_FRAMES)) {
  ref_frame <- INTERN_FRAMES[i]
  lab       <- INTERNALIZATION_LABELS[i]
  lo        <- ref_frame + 1L
  hi        <- min(ref_frame + WIN_FRAMES, MAX_FRAME)
  cat(sprintf("  intern panel %d/%d: ref t=%d (%s)  forward window frames %d-%d\n",
              i, length(INTERN_FRAMES), ref_frame, lab, lo, hi))
  fl <- build_window_flow(vel, ref_frame = ref_frame)
  if (is.null(fl)) {
    warning(sprintf("internalization panel %d (%s) had <200 spots -- skipping",
                    i, lab))
    next
  }
  fl[, anchor    := "internalization_start"]
  fl[, landmark := lab]
  flow_panels[[length(flow_panels) + 1L]] <- fl
  panel_meta[[length(panel_meta) + 1L]] <- data.table(
    anchor     = "internalization_start",
    ref_frame  = ref_frame,
    landmark   = lab,
    frame_lo   = lo,
    frame_hi   = hi,
    n_bins     = nrow(fl)
  )
}
flow_panels <- rbindlist(flow_panels, fill = TRUE)
panel_meta  <- rbindlist(panel_meta)

fwrite(flow_panels, file.path(OUT_DIR, "flow_fields.csv"))
fwrite(panel_meta,  file.path(OUT_DIR, "panel_metadata.csv"))

# Per-panel majority move_type
panel_majority <- flow_panels[, {
  n   <- .N
  maj <- names(sort(table(move_type), decreasing = TRUE))[1]
  pct <- max(table(move_type)) / n * 100
  .(majority = maj, maj_pct = pct, n_bins = n)
}, by = .(anchor, ref_frame, landmark, frame_lo, frame_hi)]
panel_majority[, t_lo_min := as.integer(frame_lo * FI_MIN)]
panel_majority[, t_hi_min := as.integer((frame_hi + 1L) * FI_MIN)]
panel_majority[, panel_lab := sprintf(
  "%s\nt=%d-%d min",
  landmark, t_lo_min, t_hi_min)]
panel_majority[, panel_lab := ascii_safe(panel_lab)]
panel_majority <- panel_majority[order(anchor, ref_frame)]
fwrite(panel_majority, file.path(OUT_DIR, "panel_majority_movement.csv"))

# Attach panel_lab and order facets by ref_frame
panel_majority[, anchor := as.character(anchor)]
flow_panels[, anchor := as.character(anchor)]
flow_panels <- merge(flow_panels,
                     panel_majority[, .(anchor, ref_frame, panel_lab)],
                     by = c("anchor", "ref_frame"))
flow_panels[, panel_lab := factor(
  panel_lab,
  levels = panel_majority[order(ref_frame), panel_lab])]

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

cat("  Plotting forward-window flow panels...\n")

plot_panel <- function(dt_sp, title_txt, subtitle_txt, ncol = 3) {
  ggplot(dt_sp) +
    geom_segment(aes(x = x_bin, y = y_bin,
                     xend = x_bin + mean_vx * ARROW_SCALE,
                     yend = y_bin + mean_vy * ARROW_SCALE,
                     color = move_type),
                 arrow = arrow(length = unit(ARROW_HEAD, "cm"),
                               type = "closed",
                               ends = "last"),
                 linewidth = ARROW_LW, alpha = 0.85) +
    facet_wrap(~ panel_lab, ncol = ncol) +
    scale_color_manual(values = flow_cols, drop = FALSE, name = NULL) +
    scale_y_reverse() + coord_fixed() +
    labs(title = title_txt, subtitle = subtitle_txt,
         x = "X (um)", y = "Y (um)  [AP top, VP bottom]") +
    theme_pub(12) +
    theme(strip.text = element_text(size = 9, face = "bold", lineheight = 1.1),
          axis.text  = element_text(size = 8),
          axis.title = element_text(size = 9))
}

n_intern <- uniqueN(flow_panels$ref_frame)

sub_short <- function() {
  sprintf("%d min AFTER each internalisation reference | %d um spatial bin | %d tracks",
          WIN_MIN, FLOW_BIN_UM, uniqueN(sp$TRACK_ID))
}

panel_intern <- plot_panel(flow_panels,
                           sprintf("mk2508 -- %d min flow AFTER internalisation landmarks", WIN_MIN),
                           sub_short(),
                           ncol = n_intern)

PANEL_W <- 5
PANEL_H <- 6

save_pdf(panel_intern,
         sprintf("04a_flow_%dmin_after_internalization_mk2508.pdf", WIN_MIN),
         w = PANEL_W * n_intern,
         h = PANEL_H)

# -----------------------------------------------------------------------------
# Combined majority-movement strip
# -----------------------------------------------------------------------------

strip_df <- copy(panel_majority)
strip_df[, anchor_label := sprintf("internalisation start (ref t=%d)",
                                   INTERNALIZATION_REFERENCE_FRAMES$internalization_start)]
strip_df[, anchor_label := factor(anchor_label, levels = unique(anchor_label))]
strip_df <- strip_df[order(anchor_label, ref_frame)]

maj_strip <- ggplot(strip_df, aes(reorder(landmark, ref_frame), majority, fill = majority)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", maj_pct)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  facet_wrap(~ anchor_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = flow_cols, drop = FALSE) +
  labs(title = sprintf("mk2508 -- majority movement per internalisation landmark, %d min forward",
                       WIN_MIN),
       subtitle = sprintf("%d-min forward window from each landmark (frames ref+1 .. ref+%d)",
                          WIN_MIN, WIN_FRAMES),
       x = "Landmark", y = NULL, fill = NULL) +
  theme_pub() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 30, 10, 10))

save_pdf(maj_strip, "04c_majority_movement_mk2508_forward.pdf",
         w = 14, h = 4)

# -----------------------------------------------------------------------------
# Audit log
# -----------------------------------------------------------------------------

ref_log_path <- file.path(OUT_DIR, "summary_references.txt")
hi_for <- function(rf) min(rf + WIN_FRAMES, MAX_FRAME)

ref_lines <- c(
  "mk2508 -- forward-window flow-map analysis (30 min AFTER each internalisation landmark)",
  "================================================================",
  sprintf("Voxel size:            %.3f um isotropic", VOXEL_UM),
  sprintf("Frame interval:        %d s (= %.2f min/frame)", FI_SEC, FI_MIN),
  sprintf("Temporal window:       %d min (= %d frames) STARTING +1 frame AFTER each landmark",
          WIN_MIN, WIN_FRAMES),
  sprintf("Velocity lag:          %d frames (= %.2f min)", VEL_LAG, VEL_LAG * FI_MIN),
  sprintf("Smoothing window:      %d frames", SMOOTH_K),
  sprintf("Spatial bin:           %d um,  min n/bin = %d",
          FLOW_BIN_UM, FLOW_MIN_N),
  sprintf("Available recording:   frames %d-%d (= %.1f min = %.2f h)",
          min(vel$FRAME), MAX_FRAME, MAX_FRAME * FI_MIN,
          MAX_FRAME * FI_MIN / 60),
  "",
  "Internalisation anchors (3 landmarks, forward 30 min):",
  sprintf("  1 h before internalisation          (ref t=%d, forward window frames %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h + 1L,
          hi_for(INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h)),
  sprintf("  Start of internalisation            (ref t=%d, forward window frames %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start + 1L,
          hi_for(INTERNALIZATION_REFERENCE_FRAMES$internalization_start)),
  sprintf("  1 h after internalisation           (ref t=%d, forward window frames %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h + 1L,
          hi_for(INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h)),
  "",
  "Outputs in this folder:",
  sprintf("  04a_flow_%dmin_after_internalization_mk2508.pdf  (3 panels)", WIN_MIN),
  "  04c_majority_movement_mk2508_forward.pdf           (majority move-type strip)",
  "  flow_fields.csv                  (raw bin-level flow data)",
  "  panel_metadata.csv               (per-panel frame window and landmark)",
  "  panel_majority_movement.csv      (per-panel majority move_type)",
  "  summary_references.txt           (this file)"
)
writeLines(ref_lines, ref_log_path)
cat("\n"); cat(paste(" ", ref_lines, collapse = "\n"), "\n", sep = "")

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
