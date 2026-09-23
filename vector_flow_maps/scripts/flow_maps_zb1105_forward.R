# =============================================================================
# FLOW MAPS — zb1105 (oriented + filtered tracks, internalisation reference)
# =============================================================================
# Mirror of the Fig-4 hourly flow panels from flow_maps_mk2508_stephane.R,
# applied to the zebrafish 1105 recording.  Zebrafish-specific changes:
#
#   * 1.25 um isotropic voxels in (z, y, x)         (vs 0.64 um for medaka)
#   * 120 s per frame                                (vs 30 s for medaka)
#   * Only the internalisation reference is used    (no contraction plot)
#   * Reference frames (absolute, supplied by the user):
#       f_1h_before_internalisation = 0
#       f_start_internalisation     = 30
#       f_1h_after_internalisation  = 60
#   * Each flow panel averages the WIN_MIN minutes STARTING AT its
#     reference frame (frame range = [ref_frame, ref_frame + WIN_FRAMES -
#     1], clipped at max_frame).  This is the opposite convention from the
#     medaka script, where the window led UP TO the reference; the flip is
#     required because the first zebrafish reference starts at frame 0, so
#     there is no data before it.
#
# Inputs (in inputs/tracks_for_flow_maps/):
#   tracks_zb_1105_oriented_filtered.csv
#
# Outputs (in outputs/zb1105/):
#   04a_flow_<WIN_MIN>min_after_internalization_zb1105.pdf   (3 panels)
#   04c_majority_movement_zb1105.pdf                         (majority strip)
#   flow_fields.csv                (binned, smoothed, vorticity-typed)
#   panel_metadata.csv             (per-panel frame window and landmark)
#   panel_majority_movement.csv    (per-panel majority move_type)
#   summary_references.txt         (audit log of reference frames used)
#
# Facet-label convention:
#   "<landmark>\nt=<lo>-<hi> min"
#   One panel per landmark.  Each panel averages every track whose frame
#   falls inside the WIN_MIN-minute window starting at the landmark frame.
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

OUT_DIR  <- "outputs/zb1105"
IN_CSV   <- file.path("inputs", "tracks_for_flow_maps",
                      "tracks_zb_1105_oriented_filtered.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Voxel calibration supplied by user: 1.25 um isotropic in (z, y, x).
VOXEL_UM <- 1.25
FI_SEC   <- 120                # frame interval (seconds/frame)
FI_MIN   <- FI_SEC / 60        # minutes per frame (= 2)

# Each flow panel averages the WIN_MIN minutes STARTING AT its reference
# frame (i.e. the WIN_MIN minutes AFTER the reference).  Frame window =
# [ref_frame, ref_frame + WIN_FRAMES - 1] inclusive, clipped at max_frame.
#
# Override via command line:  Rscript flow_maps_zb1105_stephane.R 10
# (default = 30).
WIN_MIN    <- as.integer(commandArgs(trailingOnly = TRUE)[1] %||% 30)
WIN_FRAMES <- as.integer(WIN_MIN / FI_MIN)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# Reference frames (absolute frame units) and the landmark label that goes
# with each frame.  Only the internalisation anchor is used here (the
# recording does not span contraction).
INTERNALIZATION_REFERENCE_FRAMES <- list(
  internalization_minus_1h = 0L,
  internalization_start    = 30L,
  internalization_plus_1h  = 60L
)
INTERNALIZATION_LABELS <- c(
  "1 h before internalisation",
  "Start of internalisation",
  "1 h after internalisation"
)

# Flow-field tuning (variable names match flow_maps_mk2508_stephane.R, but the
# numerical values are rescaled so the same physical interpretation holds
# at the zebrafish 120 s / frame rate: VEL_LAG = 1 frame = 2 min lag (vs
# 4 * 0.5 min = 2 min in medaka) and SMOOTH_K = 1 frame = 2 min smoothing
# (vs 5 * 0.5 min = 2.5 min in medaka).  Using medaka's frame-based values
# (4 / 5) here would give an 8 / 10 min lag/smoothing window, which is too
# long for the WIN_MIN=10 panels and would skip the frame-0 panel outright
# because the velocity computation needs more frames than the 5-frame window
# contains).
FLOW_BIN_UM       <- 30
FLOW_MIN_N        <- 10
ARROW_SCALE       <- 60
ARROW_HEAD        <- 0.14
ARROW_LW          <- 0.7
VORT_SWIRL_THRESH <- 0.3
PANEL_W <- 5
PANEL_H <- 6
SMOOTH_K          <- 1L     # 1 frame = 2 min smoothing (~medaka's 2.5 min)
VEL_LAG           <- 1L     # 1 frame = 2 min velocity lag (= medaka's 2 min)
MIN_FRAMES        <- 20L    # keep tracks with >=20 frames

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
  # The default pdf device on this R build uses a Type-1 font that can't
  # encode Unicode (Delta, minus, em-dash), so we strip non-ASCII from all
  # plot text first. The CSV outputs preserve Unicode; only the rendered PDF
  # labels are ASCII-fied.
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
# Load + orient (already oriented in the input, just apply voxel scaling)
# -----------------------------------------------------------------------------

banner("LOAD ORIENTED + FILTERED TRACKS")

sp <- fread(IN_CSV, showProgress = FALSE)
cat(sprintf("  loaded %s rows, %d columns\n", format(nrow(sp), big.mark = ","),
            ncol(sp)))
cat(sprintf("  columns: %s\n", paste(names(sp), collapse = ", ")))

# FILTERED_OUT == TRUE means "this row was dropped by the filters".
# The zebrafish oriented_filtered export in tracks_for_flow_maps/ was found
# to flag every row as FILTERED_OUT = TRUE (an export quirk), so honouring
# the flag as-is would yield an empty dataset.  When that happens we bypass
# the filter and warn; otherwise we honour it like the medaka script does.
if ("FILTERED_OUT" %in% names(sp)) {
  drop_mask <- tolower(as.character(sp$FILTERED_OUT)) == "true"
  n_in   <- nrow(sp)
  n_drop <- sum(drop_mask)
  n_keep <- n_in - n_drop
  if (n_keep == 0L && n_in > 0L) {
    warning(sprintf(
      "FILTERED_OUT flag would drop all %s rows -- bypassing filter.",
      format(n_in, big.mark = ",")))
    sp[, FILTERED_OUT := NULL]
  } else {
    cat(sprintf("  honouring FILTERED_OUT flag: %s in, %s kept, %s dropped\n",
                format(n_in, big.mark = ","),
                format(n_keep, big.mark = ","),
                format(n_drop, big.mark = ",")))
    sp <- sp[!drop_mask]
    sp[, FILTERED_OUT := NULL]
  }
}

for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z",
              "RADIAL_DIST", "SPHERICAL_DEPTH")) {
  sp[[col]] <- sp[[col]] * VOXEL_UM
}

# Short-track filter
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
cat(sprintf("  per-step rows: %s (lag = %d frames = %d min)\n",
            format(nrow(vel), big.mark = ","), VEL_LAG, VEL_LAG * FI_MIN))

# -----------------------------------------------------------------------------
# Hourly flow-field construction
# -----------------------------------------------------------------------------

banner("BUILD HOURLY FLOW FIELDS")

build_window_flow <- function(vel, ref_frame, win_frames = WIN_FRAMES,
                              bin_um = FLOW_BIN_UM, min_n = FLOW_MIN_N,
                              swirl_thr = VORT_SWIRL_THRESH,
                              max_frame = max(vel$FRAME)) {
  # Average the motion of every spot whose frame is in the window
  # [ref_frame, ref_frame + win_frames - 1] (the WIN_MIN minutes STARTING
  # AT the reference frame).  Clipped to [0, max_frame] so a window that
  # would run past the end of the recording is truncated and the panel is
  # just shorter than requested.
  frame_lo <- ref_frame
  frame_hi <- min(ref_frame + win_frames - 1L, max_frame)
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
  fl[, frame_lo := frame_lo]
  fl[, frame_hi := frame_hi]
  fl[, ref_frame := ref_frame]
  fl[]
}

# Build one flow-field per landmark.  Order matters: we walk the references
# from earliest to latest so each panel is a contiguous slice of the recording.
flow_panels <- list()
panel_meta  <- list()

INTERN_FRAMES <- as.integer(unlist(INTERNALIZATION_REFERENCE_FRAMES))

for (i in seq_along(INTERN_FRAMES)) {
  ref_frame <- INTERN_FRAMES[i]
  lab       <- INTERNALIZATION_LABELS[i]
  cat(sprintf("  intern panel %d/%d: ref t=%d (%s)  window frames %d-%d\n",
              i, length(INTERN_FRAMES), ref_frame, lab,
              ref_frame, ref_frame + WIN_FRAMES - 1L))
  fl <- build_window_flow(vel, ref_frame = ref_frame)
  if (is.null(fl)) {
    warning(sprintf("internalization panel %d (%s) had <200 spots -- skipping",
                    i, lab))
    next
  }
  fl[, landmark := lab]
  flow_panels[[length(flow_panels) + 1L]] <- fl
  panel_meta[[length(panel_meta) + 1L]] <- data.table(
    ref_frame  = ref_frame,
    landmark   = lab,
    frame_lo   = ref_frame,
    frame_hi   = min(ref_frame + WIN_FRAMES - 1L, max(vel$FRAME)),
    n_bins     = nrow(fl)
  )
}
flow_panels <- rbindlist(flow_panels, fill = TRUE)
panel_meta  <- rbindlist(panel_meta)

fwrite(flow_panels, file.path(OUT_DIR, "flow_fields.csv"))
fwrite(panel_meta,  file.path(OUT_DIR, "panel_metadata.csv"))

# Per-panel majority
panel_majority <- flow_panels[, {
  n <- .N
  maj <- names(sort(table(move_type), decreasing = TRUE))[1]
  pct <- max(table(move_type)) / n * 100
  .(majority = maj, maj_pct = pct, n_bins = n)
}, by = .(ref_frame, landmark, frame_lo, frame_hi)]
panel_majority[, t_lo_min := frame_lo * FI_MIN]
panel_majority[, t_hi_min := (frame_hi + 1L) * FI_MIN]
panel_majority[, panel_lab := sprintf(
  "%s\nt=%d-%d min",
  landmark, t_lo_min, t_hi_min)]
panel_majority[, panel_lab := ascii_safe(panel_lab)]
# Order panels by absolute reference frame so the PDF reads left-to-right
# in time order.
panel_majority <- panel_majority[order(ref_frame)]
fwrite(panel_majority, file.path(OUT_DIR, "panel_majority_movement.csv"))

# Attach the panel_lab to each row (so each facet carries its own label).
flow_panels <- merge(flow_panels,
                     panel_majority[, .(ref_frame, panel_lab)],
                     by = "ref_frame")
ord <- unique(panel_majority[order(ref_frame), panel_lab])
flow_panels[, panel_lab := factor(panel_lab, levels = ord)]

plot_panel <- function(dt_sp, title_txt, subtitle_txt, ncol = 3) {
  # Arrow length and width are pinned to PHYSICAL units (mm) so they look
  # identical regardless of how many facets share the page.
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

cat("  Plotting landmark-window flow panels...\n")

n_intern <- uniqueN(flow_panels$ref_frame)

# Subtitle: keep it to one short line of the essentials.
sub_short <- function() {
  sprintf("%d min after the reference time | %d um spatial bin | %d tracks",
          WIN_MIN, FLOW_BIN_UM, uniqueN(sp$TRACK_ID))
}

panel_intern <- plot_panel(flow_panels,
                           sprintf("zb1105 -- %d min flow vs internalisation landmarks", WIN_MIN),
                           sub_short(),
                           ncol = n_intern)

# Sizing: 3 internalisation panels in a row, 5 in x 6 in per panel.
PANEL_W <- 5
PANEL_H <- 6

save_pdf(panel_intern, sprintf("04a_flow_%dmin_after_internalization_zb1105.pdf", WIN_MIN),
         w = PANEL_W * n_intern,
         h = PANEL_H)

# -----------------------------------------------------------------------------
# Combined majority-movement strip (one row per landmark)
# -----------------------------------------------------------------------------

strip_df <- copy(panel_majority)
strip_df <- strip_df[order(ref_frame)]

maj_strip <- ggplot(strip_df, aes(reorder(landmark, ref_frame), majority, fill = majority)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", maj_pct)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = flow_cols, drop = FALSE) +
  labs(title = "zb1105 -- majority movement per landmark",
       subtitle = sprintf("%d min after the reference time", WIN_MIN),
       x = "Landmark", y = NULL, fill = NULL) +
  theme_pub() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 30, 10, 10))

save_pdf(maj_strip, "04c_majority_movement_zb1105.pdf",
         w = 14, h = 7)

# -----------------------------------------------------------------------------
# Audit log
# -----------------------------------------------------------------------------

ref_log_path <- file.path(OUT_DIR, "summary_references.txt")
ref_lines <- c(
  "zb1105 -- reference frames used in flow-map analysis",
  "================================================================",
  sprintf("Voxel size:            %.3f um isotropic", VOXEL_UM),
  sprintf("Frame interval:        %d s (= %.2f min)", FI_SEC, FI_MIN),
  sprintf("Temporal window:       %d min (= %d frames) starting at each landmark",
          WIN_MIN, WIN_FRAMES),
  sprintf("Velocity lag:          %d frames (= %.2f min)", VEL_LAG, VEL_LAG * FI_MIN),
  sprintf("Smoothing window:      %d frames (= %.2f min)", SMOOTH_K, SMOOTH_K * FI_MIN),
  sprintf("Spatial bin:           %d um,  min n/bin = %d",
          FLOW_BIN_UM, FLOW_MIN_N),
  "",
  "Anchored to START OF INTERNALIZATION (3 landmarks, absolute frames):",
  sprintf("  1 h before internalisation          (t = %d, window %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h + WIN_FRAMES - 1L),
  sprintf("  Start of internalisation            (t = %d, window %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start + WIN_FRAMES - 1L),
  sprintf("  1 h after internalisation           (t = %d, window %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h + WIN_FRAMES - 1L),
  "",
  "Outputs in this folder:",
  sprintf("  04a_flow_%dmin_after_internalization_zb1105.pdf   (3 internalisation-anchored panels)", WIN_MIN),
  "  04c_majority_movement_zb1105.pdf                   (majority move-type strip per landmark)",
  "  flow_fields.csv                (raw bin-level flow data, all panels)",
  "  panel_metadata.csv             (per-panel frame window and landmark)",
  "  panel_majority_movement.csv    (per-panel majority move_type)",
  "  summary_references.txt         (this file)"
)
writeLines(ref_lines, ref_log_path)
cat("\n"); cat(paste(" ", ref_lines, collapse = "\n"), "\n", sep = "")

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
