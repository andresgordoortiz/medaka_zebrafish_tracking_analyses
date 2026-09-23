# =============================================================================
# FLOW MAPS — mk2508 (oriented + filtered tracks)
# =============================================================================
# Port of the Fig-4 hourly flow panels from oneshot_comparison.R, applied to a
# new medaka recording (mk2508) whose tracks have already been filtered,
# sphere-fit, and oriented onto (theta, phi, depth).
#
# Reference frames supplied by the user:
#
#   Anchored to FIRST CONTRACTION (5-landmark scheme, in absolute frame units):
#     f_first_contraction       = 476   (bulge formed; recording end)
#     f_1h_before_contraction   = 356
#     f_2h_before_contraction   = 236
#     f_start_internalization   = 160
#     f_3h_before_contraction   = 116
#
#   Anchored to START OF INTERNALIZATION (3-landmark scheme):
#     f_1h_before_internalization = 40
#     f_start_internalization     = 160
#     f_1h_after_internalization  = 280
#
# Inputs (in inputs/tracks_for_flow_maps/):
#   tracks_mk_2508_stephane_oriented_filtered.csv
#
# Outputs (in outputs/mk2508/):
#   04a_flow_<WIN_MIN>min_before_landmarks_mk2508.pdf         (8 panels, all landmarks)
#   04a_flow_<WIN_MIN>min_before_contraction_mk2508.pdf       (5 panels vs first contraction)
#   04a_flow_<WIN_MIN>min_before_internalization_mk2508.pdf   (3 panels vs start of internalization)
#   04c_majority_movement_mk2508.pdf                   (majority move-type strip per landmark)
#   flow_fields.csv                (binned, smoothed, vorticity-typed)
#   panel_metadata.csv             (per-panel frame window and landmark)
#   panel_majority_movement.csv    (per-panel majority move_type)
#   summary_references.txt         (audit log of which reference frames were used)
#
# Facet-label convention:
#   "<landmark>\nt=<lo>-<hi> min"
#   One panel per landmark.  Each panel averages every track whose frame
#   falls inside the WIN_MIN-minute window ending at the landmark frame
#   (frame range = [max(0, ref_frame - WIN_FRAMES), ref_frame - 1]).
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

OUT_DIR  <- "outputs/mk2508"
IN_CSV   <- file.path("inputs", "tracks_for_flow_maps",
                      "tracks_mk_2508_stephane_oriented_filtered.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Voxel calibration supplied by user: 0.64 µm isotropic in (z, y, x).
VOXEL_UM <- 0.64
FI_SEC   <- 30                # frame interval (seconds/frame)
FI_MIN   <- FI_SEC / 60       # minutes per frame

# Each flow panel averages the WIN_MIN minutes LEADING UP TO its reference
# frame (i.e. the WIN_MIN minutes ending at the landmark).  Frame window =
# [t - WIN, t - 1] inclusive so the panel ends exactly at the landmark frame.
#
# Override via command line:  Rscript flow_maps_mk2508_stephane.R 10
# (default = 30).
WIN_MIN <- as.integer(commandArgs(trailingOnly = TRUE)[1] %||% 30)
WIN_FRAMES <- as.integer(WIN_MIN / FI_MIN)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# Reference frames (absolute frame units) and the landmark label that goes
# with each frame.  The contraction panel reports the closest landmark to the
# bin midpoint (5 landmarks); the internalization panel reports the closest
# landmark (3 landmarks).
CONTRACTION_REFERENCE_FRAMES <- list(
  contraction_first      = 476L,
  contraction_minus_1h   = 356L,
  contraction_minus_2h   = 236L,
  internalization_start  = 160L,
  contraction_minus_3h   = 116L
)
CONTRACTION_LABELS <- c(
  "First contraction / bulge formed",
  "1 h before first contraction",
  "2 h before first contraction",
  "Start of internalization",
  "3 h before first contraction"
)
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

# Flow-field tuning
FLOW_BIN_UM       <- 30
FLOW_MIN_N        <- 10
ARROW_SCALE       <- 60
ARROW_HEAD        <- 0.14
ARROW_LW          <- 0.7
VORT_SWIRL_THRESH <- 0.3
PANEL_W <- 5
PANEL_H <- 6
SMOOTH_K          <- 5L
VEL_LAG           <- 4L     # 2-minute step
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
                              swirl_thr = VORT_SWIRL_THRESH) {
  # Average the motion of every spot whose frame is in the window
  # [ref_frame - win_frames, ref_frame - 1] (the WIN_MIN minutes leading up
  # to the reference frame).  If the window starts before frame 0 (i.e. the
  # reference is too early in the recording), clip to [0, ref_frame - 1] and
  # return the effective window in frame_lo / frame_hi so the caller knows
  # the panel is shorter than requested.
  frame_hi <- ref_frame - 1L
  frame_lo <- max(0L, ref_frame - win_frames)
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

CONTR_FRAMES  <- as.integer(unlist(CONTRACTION_REFERENCE_FRAMES))
INTERN_FRAMES <- as.integer(unlist(INTERNALIZATION_REFERENCE_FRAMES))

for (i in seq_along(CONTR_FRAMES)) {
  ref_frame <- CONTR_FRAMES[i]
  lab       <- CONTRACTION_LABELS[i]
  cat(sprintf("  contr panel %d/%d: ref t=%d (%s)  window frames %d-%d\n",
              i, length(CONTR_FRAMES), ref_frame, lab,
              ref_frame - WIN_FRAMES, ref_frame - 1L))
  fl <- build_window_flow(vel, ref_frame = ref_frame)
  if (is.null(fl)) {
    warning(sprintf("contraction panel %d (%s) had <200 spots -- skipping",
                    i, lab))
    next
  }
  fl[, anchor    := "first_contraction"]
  fl[, landmark := lab]
  flow_panels[[length(flow_panels) + 1L]] <- fl
  panel_meta[[length(panel_meta) + 1L]] <- data.table(
    anchor     = "first_contraction",
    ref_frame  = ref_frame,
    landmark   = lab,
    frame_lo   = ref_frame - WIN_FRAMES,
    frame_hi   = ref_frame - 1L,
    n_bins     = nrow(fl)
  )
}
for (i in seq_along(INTERN_FRAMES)) {
  ref_frame <- INTERN_FRAMES[i]
  lab       <- INTERNALIZATION_LABELS[i]
  cat(sprintf("  intern panel %d/%d: ref t=%d (%s)  window frames %d-%d\n",
              i, length(INTERN_FRAMES), ref_frame, lab,
              ref_frame - WIN_FRAMES, ref_frame - 1L))
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
    frame_lo   = ref_frame - WIN_FRAMES,
    frame_hi   = ref_frame - 1L,
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
}, by = .(anchor, ref_frame, landmark, frame_lo, frame_hi)]
panel_majority[, t_lo_min := frame_lo * FI_MIN]
panel_majority[, t_hi_min := (frame_hi + 1L) * FI_MIN]
# Compact facet header: "<landmark> | t<frame_lo>-<frame_hi>"
# ("maj:%s %.0f%%" is dropped from the strip text -- it's already on the
# bottom strip chart.)
panel_majority[, panel_lab := sprintf(
  "%s\nt=%d-%d min",
  landmark, t_lo_min, t_hi_min)]
panel_majority[, panel_lab := ascii_safe(panel_lab)]
# Order panels by absolute reference frame so each PDF reads left-to-right
# in time order.
panel_majority <- panel_majority[order(anchor, ref_frame)]
fwrite(panel_majority, file.path(OUT_DIR, "panel_majority_movement.csv"))

# Attach the panel_lab to each row (so each facet carries its own label).
panel_majority[, anchor := as.character(anchor)]
flow_panels[, anchor := as.character(anchor)]
flow_panels <- merge(flow_panels,
                     panel_majority[, .(anchor, ref_frame, panel_lab)],
                     by = c("anchor", "ref_frame"))

# Split by anchor for the two PDFs, and assign each subset's panel_lab
# levels in absolute-frame order so the facets render left-to-right in time.
set_panel_order <- function(dt, anchor_name) {
  sub <- dt[anchor == anchor_name]
  ord <- unique(panel_majority[anchor == anchor_name][order(ref_frame), panel_lab])
  sub[, panel_lab := factor(panel_lab, levels = ord)]
  sub[]
}
flow_contr  <- set_panel_order(flow_panels, "first_contraction")
flow_intern <- set_panel_order(flow_panels, "internalization_start")
# All-landmarks view: both anchors, ordered by absolute ref_frame.
flow_all <- copy(flow_panels)
all_ord <- unique(panel_majority[order(ref_frame), panel_lab])
flow_all[, panel_lab := factor(panel_lab, levels = all_ord)]

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

n_contr  <- uniqueN(flow_contr$ref_frame)
n_intern <- uniqueN(flow_intern$ref_frame)
n_all    <- uniqueN(flow_all$ref_frame)

# Subtitle: keep it to one short line of the essentials.
sub_short <- function() {
  sprintf("%d min before the reference time | %d um spatial bin | %d tracks",
          WIN_MIN, FLOW_BIN_UM, uniqueN(sp$TRACK_ID))
}

panel_contr <- plot_panel(flow_contr,
                          sprintf("mk2508 -- %d min flow vs first-contraction landmarks", WIN_MIN),
                          sub_short(),
                          ncol = n_contr)   # one row of 5 columns

panel_intern <- plot_panel(flow_intern,
                           sprintf("mk2508 -- %d min flow vs internalization landmarks", WIN_MIN),
                           sub_short(),
                           ncol = n_intern)

panel_all <- plot_panel(flow_all,
                        sprintf("mk2508 -- %d min flow, all landmarks", WIN_MIN),
                        sub_short(),
                        ncol = if (n_all <= 4) n_all else 4)

# Sizing: keep the internalization PDF exactly as it was (3 panels in a
# row, 15 in x 6 in) and match the contraction PDF to that exact panel
# size: 5 in wide, 6 in tall per panel.
PANEL_W <- 5
PANEL_H <- 6

save_pdf(panel_contr,  sprintf("04a_flow_%dmin_before_contraction_mk2508.pdf", WIN_MIN),
         w = PANEL_W * n_contr,
         h = PANEL_H)
save_pdf(panel_intern, sprintf("04a_flow_%dmin_before_internalization_mk2508.pdf", WIN_MIN),
         w = PANEL_W * n_intern,
         h = PANEL_H)
save_pdf(panel_all,    sprintf("04a_flow_%dmin_before_landmarks_mk2508.pdf", WIN_MIN),
         w = PANEL_W * if (n_all <= 4) n_all else 4,
         h = PANEL_H * ceiling(n_all    / if (n_all    <= 4) n_all    else 4))

# -----------------------------------------------------------------------------
# Combined majority-movement strip (one row per landmark, faceted by anchor)
# -----------------------------------------------------------------------------

strip_df <- copy(panel_majority)
strip_df[, anchor_label := fifelse(anchor == "first_contraction",
                                   sprintf("first contraction (ref t=%d)",
                                           CONTRACTION_REFERENCE_FRAMES$contraction_first),
                                   sprintf("internalization start (ref t=%d)",
                                           INTERNALIZATION_REFERENCE_FRAMES$internalization_start))]
strip_df[, anchor_label := factor(anchor_label, levels = unique(anchor_label))]
strip_df <- strip_df[order(anchor_label, ref_frame)]

maj_strip <- ggplot(strip_df, aes(reorder(landmark, ref_frame), majority, fill = majority)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", maj_pct)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  facet_wrap(~ anchor_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = flow_cols, drop = FALSE) +
  labs(title = "mk2508 -- majority movement per landmark",
       subtitle = sprintf("%d min before the reference time", WIN_MIN),
       x = "Landmark", y = NULL, fill = NULL) +
  theme_pub() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 30, 10, 10))

save_pdf(maj_strip, "04c_majority_movement_mk2508.pdf",
         w = 14, h = 7)

# -----------------------------------------------------------------------------
# Audit log
# -----------------------------------------------------------------------------

ref_log_path <- file.path(OUT_DIR, "summary_references.txt")
ref_lines <- c(
  "mk2508 -- reference frames used in flow-map analysis",
  "================================================================",
  sprintf("Voxel size:            %.3f um isotropic", VOXEL_UM),
  sprintf("Frame interval:        %d s (= %.2f min)", FI_SEC, FI_MIN),
  sprintf("Temporal window:       %d min (= %d frames) ending at each landmark",
          WIN_MIN, WIN_FRAMES),
  sprintf("Velocity lag:          %d frames (= %.2f min)", VEL_LAG, VEL_LAG * FI_MIN),
  sprintf("Smoothing window:      %d frames", SMOOTH_K),
  sprintf("Spatial bin:           %d um,  min n/bin = %d",
          FLOW_BIN_UM, FLOW_MIN_N),
  "",
  "Anchored to FIRST CONTRACTION (5 landmarks, absolute frames):",
  sprintf("  3 h before first contraction        (t = %d, window %d-%d)",
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_3h,
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_3h - WIN_FRAMES,
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_3h - 1L),
  sprintf("  Start of internalization            (t = %d, window %d-%d)",
          CONTRACTION_REFERENCE_FRAMES$internalization_start,
          CONTRACTION_REFERENCE_FRAMES$internalization_start - WIN_FRAMES,
          CONTRACTION_REFERENCE_FRAMES$internalization_start - 1L),
  sprintf("  2 h before first contraction        (t = %d, window %d-%d)",
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_2h,
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_2h - WIN_FRAMES,
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_2h - 1L),
  sprintf("  1 h before first contraction        (t = %d, window %d-%d)",
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_1h,
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_1h - WIN_FRAMES,
          CONTRACTION_REFERENCE_FRAMES$contraction_minus_1h - 1L),
  sprintf("  First contraction / bulge formed    (t = %d, window %d-%d)",
          CONTRACTION_REFERENCE_FRAMES$contraction_first,
          CONTRACTION_REFERENCE_FRAMES$contraction_first - WIN_FRAMES,
          CONTRACTION_REFERENCE_FRAMES$contraction_first - 1L),
  "",
  "Anchored to START OF INTERNALIZATION (3 landmarks, absolute frames):",
  sprintf("  1 h before internalization          (t = %d, window %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h - WIN_FRAMES,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_minus_1h - 1L),
  sprintf("  Start of internalization            (t = %d, window %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start - WIN_FRAMES,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_start - 1L),
  sprintf("  1 h after internalization           (t = %d, window %d-%d)",
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h - WIN_FRAMES,
          INTERNALIZATION_REFERENCE_FRAMES$internalization_plus_1h - 1L),
  "",
  "Outputs in this folder:",
  sprintf("  04a_flow_%dmin_before_landmarks_mk2508.pdf         (all 8 panels)", WIN_MIN),
  sprintf("  04a_flow_%dmin_before_contraction_mk2508.pdf       (5 contraction-anchored panels)", WIN_MIN),
  sprintf("  04a_flow_%dmin_before_internalization_mk2508.pdf   (3 internalization-anchored panels)", WIN_MIN),
  "  04c_majority_movement_mk2508.pdf                   (majority move-type strip per landmark)",
  "  flow_fields.csv                (raw bin-level flow data, all panels)",
  "  panel_metadata.csv             (per-panel frame window and landmark)",
  "  panel_majority_movement.csv    (per-panel majority move_type)",
  "  summary_references.txt         (this file)"
)
writeLines(ref_lines, ref_log_path)
cat("\n"); cat(paste(" ", ref_lines, collapse = "\n"), "\n", sep = "")

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
