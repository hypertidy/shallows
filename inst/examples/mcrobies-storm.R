## Rain on grid: a 30 minute, 60 mm/h storm on a 1.2 km tile of the
## Tasmania 2 m DEM around McRobies Gully (Hobart), then 10 minutes of
## drainage. Open boundaries on all four edges. Manning n = 0.05 everywhere
## (no land cover, no infiltration, no drains): this is a test of the
## solver on real terrain, for illustration only, not survey, engineering
## or planning data.
##
## Raster input and output with gdalraster, plots with ximage.
## Usage: Rscript mcrobies-storm.R [outdir] [aggregate_factor]
suppressMessages(library(shallows))
library(gdalraster)
library(ximage)
args <- commandArgs(TRUE)
out <- if (length(args) >= 1) args[1] else "mcrobies-storm"
fact <- if (length(args) >= 2) as.integer(args[2]) else 2L
dir.create(out, showWarnings = FALSE, recursive = TRUE)

dsn <- "/vsicurl/https://data.source.coop/alexgleith/tasmania-dem-2m/Tasmania_Statewide_2m_DEM_14-08-2021.tif"
centre <- c(523427, 5251627)  # EPSG:28355, McRobies Gully tip floor (147.2869 E, 42.8876 S)
te <- c(centre - 600, centre + 600)  # xmin, ymin, xmax, ymax
bed_file <- file.path(out, "bed.tif")
## crop, and average to fact * 2 m cells, in one gdalwarp
warp(dsn, bed_file, t_srs = "EPSG:28355", quiet = TRUE,
     cl_arg = c("-te", te, "-tr", 2 * fact, 2 * fact, "-r", "average", "-overwrite"))

rain <- 60 / 1000 / 3600  # 60 mm/h in m/s
st <- sw_from_rast(bed_file, h = 0, manning = 0.05, rain = rain)
t1 <- system.time(r1 <- sw_run(st, times = seq(60, 1800, by = 60), order = 2, bc = "open"))
st2 <- sw_continue(r1, rain = 0)
t2 <- system.time(r2 <- sw_run(st2, times = seq(60, 600, by = 60), order = 2, bc = "open"))

## depth frames with CRS, one band per minute
sw_write(r1, file.path(out, "depth_rain.tif"), "h")
sw_write(r2, file.path(out, "depth_drain.tif"), "h")
sw_write(r2, file.path(out, "speed_drain.tif"), "speed")
h1 <- sw_frames(r1, "h"); h2 <- sw_frames(r2, "h")
hmax <- Reduce(pmax, lapply(c(h1, h2), sw_matrix))

d <- rbind(transform(r1$diag, phase = "rain"),
           transform(r2$diag, phase = "drain", time = time + 1800,
                     outflow = outflow + r1$diag$outflow[30], rain = rain + r1$diag$rain[30]))
write.csv(d, file.path(out, "budget.csv"), row.names = FALSE)
cat(sprintf("cells %d x %d at %g m; steps %d + %d; wall time %.0f + %.0f s\n",
            nrow(st$z), ncol(st$z), st$dx, max(r1$diag$steps), max(r2$diag$steps), t1[3], t2[3]))
cat(sprintf("rain in %.0f m3, out through edges %.0f m3, stored %.0f m3, max relative budget error %.1e\n",
            max(d$rain), max(d$outflow), d$volume[nrow(d)],
            max(abs(c(r1$diag$rel_budget_error, r2$diag$rel_budget_error)))))

## figure: depth over a hillshade
hs_file <- file.path(out, "hillshade.tif")
dem_proc("hillshade", bed_file, hs_file, mode_options = c("-az", "300", "-alt", "40"), quiet = TRUE)
hs <- read_ds(new(GDALRaster, hs_file))
ext <- st$extent
caveat <- "For illustration only: not survey, engineering or planning data. Bed: Tasmania 2 m DEM (2021). Uniform n = 0.05, no infiltration or drains."
pal <- hcl.colors(96, "Blues 3", rev = TRUE)
show <- function(h, main) {
  ximage(hs, col = grey.colors(96, 0, 1), asp = 1, axes = FALSE, main = main)
  h[h < 0.01] <- NA
  ximage(h, extent = ext, zlim = c(0.01, 1), col = pal, alpha = 0.85, add = TRUE)
}
png(file.path(out, "mcrobies-storm.png"), width = 1800, height = 1100, res = 150)
layout(matrix(c(1, 2, 3, 4, 4, 4), 2, byrow = TRUE), heights = c(2.3, 1.5))
par(mar = c(1, 1, 3, 1))
show(sw_matrix(h1[[10]]), "depth at 10 min")
show(sw_matrix(h1[[30]]), "depth at 30 min (rain stops)")
show(hmax, "maximum depth, 0 to 40 min (0.01 to 1 m)")
par(mar = c(7, 5, 2, 5))
tt <- c(0, d$time) / 60
plot(tt, c(0, d$rain), type = "l", lty = 2, xlab = "minutes", ylab = "m3",
     main = "water budget of the tile", ylim = c(0, max(d$rain)))
lines(tt, c(0, d$volume), col = "steelblue", lwd = 2)
lines(tt, c(0, d$outflow), col = "firebrick", lwd = 2)
legend("topleft", c("rain in (cumulative)", "stored on the tile", "out through edges (cumulative)"),
       lty = c(2, 1, 1), col = c("black", "steelblue", "firebrick"), lwd = c(1, 2, 2), bty = "n")
mtext(caveat, side = 1, line = 5.2, cex = 0.75)
invisible(dev.off())
