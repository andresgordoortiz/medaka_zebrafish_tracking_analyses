suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

OUT  <- "results/oneshot_comparison"
dir.create(OUT, showWarnings = FALSE)

MED_DIR <- "data/oriented_medaka_ultrack"
MED_VOX <- 1.05152
MED_FI  <- 30  # s/frame

sp <- fread(file.path(MED_DIR, "oriented_tracks_medaka.csv"), showProgress = FALSE)
for (col in c("POSITION_X","POSITION_Y","POSITION_Z","RADIAL_DIST","SPHERICAL_DEPTH"))
  sp[[col]] <- sp[[col]] * MED_VOX
sp[, time_min := FRAME * MED_FI / 60]
R <- as.numeric(fread(file.path(MED_DIR, "sphere_params.csv"))[parameter=="radius", value]) * MED_VOX
cat(sprintf("Medaka: sphere R=%.1f um, %s spots, %s tracks, %.0f min\n",
            R, format(nrow(sp), big.mark=","), format(uniqueN(sp$TRACK_ID), big.mark=","),
            max(sp$time_min)))

# === Re-derive the bulge the same way as oneshot ============================
late <- sp[FRAME >= max(FRAME)/2 & is.finite(SPHERICAL_DEPTH)]
thr  <- quantile(late$SPHERICAL_DEPTH, 0.95)
cand <- late[SPHERICAL_DEPTH >= thr]
ct <- mean(cand$THETA_DEG); cp <- mean(cand$PHI_DEG)
cand[, ang := sqrt((THETA_DEG-ct)^2 + (PHI_DEG-cp)^2)]
rad <- min(2*sd(cand$ang)+2, 15)
cat(sprintf("\nDATA-DRIVEN BULGE:  theta=%.2f  phi=%.2f  radius=%.2f deg  (depth>=%.1f um in late half)\n",
            ct, cp, rad, thr))

# Try a narrower, depth-weighted centre: re-center on top 1% only
top1 <- sp[FRAME >= max(FRAME)/2 & SPHERICAL_DEPTH >= quantile(sp$SPHERICAL_DEPTH, 0.99, na.rm=TRUE)]
cat(sprintf("Top-1%% centre:    theta=%.2f  phi=%.2f  (n=%d cells)\n",
            mean(top1$THETA_DEG, na.rm=TRUE), mean(top1$PHI_DEG, na.rm=TRUE), nrow(top1)))

# === 1.  Where IS the maximum density on the embryo at late time? ==========
late30 <- sp[time_min >= max(time_min) - 30 & is.finite(SPHERICAL_DEPTH)]
late30[, theta_bin := round(THETA_DEG / 4) * 4]
late30[, phi_bin   := round(PHI_DEG   / 4) * 4]
heat <- late30[, .(n = .N,
                   p50_d = median(SPHERICAL_DEPTH),
                   p90_d = quantile(SPHERICAL_DEPTH, 0.90),
                   nucl_density_2D = .N),
               by = .(theta_bin, phi_bin)][n >= 5]

cat("\nTOP 10 hottest 4x4 deg bins in late 30 min (by # nuclei in 2D theta-phi projection):\n")
print(heat[order(-n)][1:10, .(theta_bin, phi_bin, n_nuclei = n, p50_depth = round(p50_d,1), p90_depth = round(p90_d,1))])

cat("\nTOP 10 DEEPEST 4x4 deg bins in late 30 min (by p90 depth):\n")
print(heat[order(-p90_d)][1:10, .(theta_bin, phi_bin, n_nuclei = n, p50_depth = round(p50_d,1), p90_depth = round(p90_d,1))])

# === 2.  Count of cells inside the disc OVER TIME, all depths and split ====
disc_cells <- function(sp, ct, cp, r) {
  sp[!is.na(PHI_DEG) & sqrt((THETA_DEG-ct)^2 + (PHI_DEG-cp)^2) < r]
}
sp_disc <- disc_cells(sp, ct, cp, rad)
cat(sprintf("\nDisc contains %s total observations (%s tracks)\n",
            format(nrow(sp_disc), big.mark=","), format(uniqueN(sp_disc$TRACK_ID), big.mark=",")))

per_frame <- sp_disc[, .(
    n_total   = .N,
    n_deep    = sum(SPHERICAL_DEPTH >= 30, na.rm = TRUE),
    n_surface = sum(SPHERICAL_DEPTH <  30, na.rm = TRUE),
    p50_depth = median(SPHERICAL_DEPTH, na.rm = TRUE),
    p90_depth = quantile(SPHERICAL_DEPTH, 0.90, na.rm = TRUE)),
  by = FRAME][order(FRAME)]
per_frame[, time_min := FRAME * MED_FI / 60]

cat("\nMedaka bulge disc -- counts vs time (snapshot every 50 frames):\n")
print(per_frame[seq(1, nrow(per_frame), by = 50)])

# Disc area
A_disc <- 2*pi*R^2*(1 - cos(rad*pi/180))
cat(sprintf("\nDisc area on sphere: %.1f x 10^3 um^2\n", A_disc/1000))

per_frame[, dens_all  := n_total   / A_disc * 1000]
per_frame[, dens_deep := n_deep    / A_disc * 1000]
per_frame[, dens_surf := n_surface / A_disc * 1000]

# Smooth and compare
plot_dt <- melt(per_frame, id.vars=c("FRAME","time_min"),
                measure.vars=c("dens_all","dens_deep","dens_surf"),
                variable.name="layer", value.name="density")
plot_dt[, layer := factor(layer, levels=c("dens_all","dens_surf","dens_deep"),
                          labels=c("All depths","Surface (depth<30um)","Deep (depth>=30um)"))]

p_dens <- ggplot(plot_dt, aes(time_min, density, color=layer)) +
  geom_point(alpha=0.25, size=0.5) +
  geom_smooth(method="loess", span=0.3, se=FALSE, linewidth=1) +
  scale_color_manual(values=c("All depths"="black","Surface (depth<30um)"="#1F78B4","Deep (depth>=30um)"="#E31A1C")) +
  labs(title=sprintf("Medaka bulge disc (theta=%.1f phi=%.1f r=%.1f deg, A=%.1f x 10^3 um^2)",
                     ct, cp, rad, A_disc/1000),
       subtitle="Counts of nuclei in the disc, by depth layer, normalised by disc area",
       x="time (min)", y="nuclei / 1000 um^2", color=NULL) +
  theme_minimal()

# === 3.  Compare disc vs a CONTROL disc on the other side ==================
ctrl_ct <- ct
ctrl_cp <- ((cp + 180) %% 360) - 180   # antipodal in phi, same theta
ctrl_disc <- disc_cells(sp, ctrl_ct, ctrl_cp, rad)
cat(sprintf("\nCONTROL disc (antipodal phi=%.1f, same theta=%.1f): %d obs\n",
            ctrl_cp, ctrl_ct, nrow(ctrl_disc)))

per_frame_c <- ctrl_disc[, .(n_total=.N, p50_depth=median(SPHERICAL_DEPTH, na.rm=TRUE)),
                          by=FRAME][order(FRAME)]
per_frame_c[, time_min := FRAME * MED_FI / 60]
per_frame_c[, dens := n_total / A_disc * 1000]

cmp <- rbind(
  data.table(time_min=per_frame$time_min,   density=per_frame$dens_all, location="Bulge"),
  data.table(time_min=per_frame_c$time_min, density=per_frame_c$dens,   location="Control (antipodal phi)"))
p_cmp <- ggplot(cmp, aes(time_min, density, color=location)) +
  geom_point(alpha=0.25, size=0.5) +
  geom_smooth(method="loess", span=0.3, se=FALSE, linewidth=1) +
  scale_color_manual(values=c("Bulge"="#E31A1C","Control (antipodal phi)"="#1F78B4")) +
  labs(title="Disc density: bulge vs control disc (same theta, opposite phi)",
       subtitle="If bulge truly thickens, RED should rise above BLUE; if both flat-and-similar, density is not picking it up",
       x="time (min)", y="nuclei / 1000 um^2", color=NULL) +
  theme_minimal()

# === 4.  Volumetric density inside vs outside the disc ====================
# Same surface area but include DEPTH -- volume = A_disc * thickness band
# For each frame compute n_in_disc and the "column volume" up to depth p99
sp_finite <- sp[is.finite(SPHERICAL_DEPTH)]
sp_finite[, in_disc := sqrt((THETA_DEG-ct)^2 + (PHI_DEG-cp)^2) < rad]
vol_dt <- sp_finite[, .(
    n_in     = sum(in_disc),
    n_out    = sum(!in_disc),
    d_in_p10 = quantile(SPHERICAL_DEPTH[in_disc],  0.10, na.rm=TRUE),
    d_in_p90 = quantile(SPHERICAL_DEPTH[in_disc],  0.90, na.rm=TRUE),
    d_out_p10= quantile(SPHERICAL_DEPTH[!in_disc], 0.10, na.rm=TRUE),
    d_out_p90= quantile(SPHERICAL_DEPTH[!in_disc], 0.90, na.rm=TRUE)),
  by=FRAME][order(FRAME)]
vol_dt[, time_min := FRAME * MED_FI / 60]
vol_dt[, V_in  := A_disc * pmax(d_in_p90  - d_in_p10,  1)]   # um^3
# global area: use crude sphere patch
A_out <- 2*pi*R^2*0.5   # roughly half sphere
vol_dt[, V_out := A_out * pmax(d_out_p90 - d_out_p10, 1)]
vol_dt[, vol_dens_in  := n_in  / V_in  * 1e6]   # nuclei per 10^6 um^3
vol_dt[, vol_dens_out := n_out / V_out * 1e6]

p_vol <- ggplot(melt(vol_dt[, .(time_min, vol_dens_in, vol_dens_out)],
                     id.vars="time_min", variable.name="loc", value.name="dens"),
                aes(time_min, dens, color=loc)) +
  geom_point(alpha=0.3, size=0.5) +
  geom_smooth(method="loess", span=0.3, se=FALSE, linewidth=1) +
  scale_color_manual(values=c("vol_dens_in"="#E31A1C","vol_dens_out"="#1F78B4"),
                     labels=c("Inside bulge disc (vol)","Outside bulge disc (vol)")) +
  labs(title="VOLUMETRIC nuclei density inside vs outside the disc (per 10^6 um^3)",
       subtitle="Using p10-p90 depth range as the column thickness (so this approximates true cells/volume)",
       x="time (min)", y="nuclei / 10^6 um^3", color=NULL) +
  theme_minimal()

# === 5.  Spatial maps at start vs end with the disc overlay ===============
make_circle <- function(ct, cp, r, n=80){
  a <- seq(0, 2*pi, length.out=n)
  data.table(phi = cp + r*cos(a), theta = ct + r*sin(a))
}
circ <- make_circle(ct, cp, rad)

snap_t <- c(round(max(sp$time_min)*0.05),
            round(max(sp$time_min)*0.5),
            round(max(sp$time_min)*0.95))
snaps <- rbindlist(lapply(snap_t, function(t)
  cbind(sp[time_min >= t-5 & time_min <= t+5, .(PHI_DEG, THETA_DEG, SPHERICAL_DEPTH)],
        snap = sprintf("t = %d min", t))))
snaps[, snap := factor(snap, levels=sprintf("t = %d min", snap_t))]

p_map <- ggplot(snaps, aes(PHI_DEG, THETA_DEG)) +
  geom_point(aes(color=SPHERICAL_DEPTH), size=0.4, alpha=0.5) +
  geom_polygon(data=circ, aes(phi, theta), fill=NA, color="red", linewidth=0.7, inherit.aes=FALSE) +
  facet_wrap(~snap, nrow=1) +
  scale_y_reverse() +
  scale_color_viridis_c(option="inferno", name="depth (um)") +
  labs(title="Medaka — spatial distribution of nuclei (sphere surface, theta vs phi)",
       subtitle="Red circle = data-driven bulge disc.  Colour = spherical depth.  Each dot = one nucleus observation.",
       x="phi (deg)", y="theta (deg)") +
  theme_minimal()

fig <- (p_map / p_dens / p_cmp / p_vol) + plot_layout(heights = c(1.5, 1, 1, 1)) +
  plot_annotation(title = "Medaka bulge density — diagnostic",
                  theme = theme(plot.title = element_text(face="bold", size=14)))
ggsave(file.path(OUT, "medaka_bulge_density_check.pdf"), fig, w=13, h=18)
cat(sprintf("\nSaved %s\n", file.path(OUT, "medaka_bulge_density_check.pdf")))

# === 6.  Quick numbers ============================
cat("\n=== Summary numbers for medaka bulge disc ===\n")
early <- per_frame[time_min <= 30]
late  <- per_frame[time_min >= max(time_min)-30]
cat(sprintf("Surface density (all depths):  early mean = %.3f, late mean = %.3f  (ratio %.2fx)\n",
            mean(early$dens_all),  mean(late$dens_all),  mean(late$dens_all)/mean(early$dens_all)))
cat(sprintf("Deep-cell density (>=30um):    early mean = %.3f, late mean = %.3f  (ratio %.2fx)\n",
            mean(early$dens_deep), mean(late$dens_deep), mean(late$dens_deep)/mean(early$dens_deep)))
cat(sprintf("Median depth in disc:          early = %.1f um, late = %.1f um\n",
            mean(early$p50_depth), mean(late$p50_depth)))
cat(sprintf("p90 depth in disc:             early = %.1f um, late = %.1f um\n",
            mean(early$p90_depth), mean(late$p90_depth)))

cat("\nVOLUMETRIC numbers (inside disc):\n")
vol_early <- vol_dt[time_min<=30]; vol_late <- vol_dt[time_min>=max(time_min)-30]
cat(sprintf("Volumetric density inside disc: early = %.3f, late = %.3f (ratio %.2fx) per 10^6 um^3\n",
            mean(vol_early$vol_dens_in), mean(vol_late$vol_dens_in),
            mean(vol_late$vol_dens_in)/mean(vol_early$vol_dens_in)))
cat(sprintf("Volumetric density OUTSIDE disc: early = %.3f, late = %.3f (ratio %.2fx) per 10^6 um^3\n",
            mean(vol_early$vol_dens_out), mean(vol_late$vol_dens_out),
            mean(vol_late$vol_dens_out)/mean(vol_early$vol_dens_out)))
