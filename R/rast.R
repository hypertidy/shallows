## Grids in and out.
##
## A grid is a numeric vector in the layout of gdalraster::read_ds(): cell
## values row major, first row at the top (north), with an attribute "gis"
## holding bbox (xmin, ymin, xmax, ymax), dim (ncol, nrow, nbands), srs (WKT,
## or "" for local metres) and datatype. ximage::ximage() plots a
## single-band grid directly, in place, with its extent.

.need_gdalraster <- function() {
  if (!requireNamespace("gdalraster", quietly = TRUE)) stop("needs the gdalraster package")
}

.is_gis <- function(x) is.list(attr(x, "gis"))

.gis <- function(m, extent, srs = "", ...) {
  m <- as.matrix(m)
  structure(as.vector(t(m)),
            gis = list(type = "raster", bbox = extent[c(1, 3, 2, 4)],
                       dim = c(ncol(m), nrow(m), 1), srs = if (is.null(srs)) "" else srs,
                       datatype = "Float64"),
            ...)
}

.gis_extent <- function(x) attr(x, "gis")$bbox[c(1, 3, 2, 4)]

.as_wkt <- function(srs) {
  if (is.null(srs) || !nzchar(srs)) return("")
  gdalraster::srs_to_wkt(srs)
}

#' Matrix from a grid
#'
#' A grid as returned by [sw_generate()], [sw_extract()], [sw_frames()] or
#' `gdalraster::read_ds()` is a row-major vector with a `gis` attribute; this
#' gives one band of it as a matrix oriented like the raster (row 1 at the
#' top), which is the orientation [sw_state()] and `ximage::ximage()` use.
#' Anything else is passed to [as.matrix()].
#' @param x grid, or anything [as.matrix()] accepts
#' @param band band number
#' @return matrix
#' @export
sw_matrix <- function(x, band = 1L) {
  if (!.is_gis(x)) return(as.matrix(x))
  g <- attr(x, "gis")
  nx <- g$dim[1]; ny <- g$dim[2]; n <- nx * ny
  matrix(as.vector(x)[(band - 1) * n + seq_len(n)], ny, nx, byrow = TRUE)
}

.read_gis <- function(x) {
  .need_gdalraster()
  if (is.character(x)) {
    x <- methods::new(gdalraster::GDALRaster, x)
    on.exit(x$close())
  }
  gdalraster::read_ds(x, bands = 1L)
}

#' Build a state from a DEM
#'
#' Reads the bed from a raster, keeping its extent and CRS so outputs can be
#' written back as georeferenced rasters. Missing bed cells are filled with
#' the maximum elevation plus 100 m (a high wall).
#' @param dem the bed in metres, in a projected CRS (or none, for a generated
#'   bed): a file or GDAL source name, a `gdalraster::GDALRaster`, or a grid
#'   as from [sw_generate()], [sw_extract()] or `gdalraster::read_ds()`. Only
#'   the first band is used
#' @param ... passed to [sw_state()] (h, eta, manning, rain, ...)
#' @export
sw_from_rast <- function(dem, ...) {
  if (is.character(dem) || inherits(dem, "Rcpp_GDALRaster")) dem <- .read_gis(dem)
  if (!.is_gis(dem))
    stop("'dem' must be a file or GDAL source name, a GDALRaster, or a grid as from read_ds() or sw_generate()")
  g <- attr(dem, "gis")
  if (nzchar(g$srs)) {
    .need_gdalraster()
    if (gdalraster::srs_is_geographic(g$srs))
      stop("DEM must be in a projected CRS with metre units (see sw_extract())")
  }
  z <- sw_matrix(dem)
  if (anyNA(z)) z[is.na(z)] <- max(z, na.rm = TRUE) + 100
  ext <- .gis_extent(dem)
  sw_state(z, dx = (ext[2] - ext[1]) / g$dim[1], dy = (ext[4] - ext[3]) / g$dim[2],
           extent = ext, crs = g$srs, ...)
}

.run_extent <- function(s) {
  if (is.null(s$extent)) c(0, ncol(s$z) * s$dx, 0, nrow(s$z) * s$dy) else s$extent
}

#' Frames of a run as georeferenced grids
#'
#' One grid per output time, each a row-major vector with a `gis` attribute
#' (extent and CRS of the state) and a `time` attribute, as read by
#' `gdalraster::read_ds()`; plot one with `ximage::ximage()`, get a matrix
#' with [sw_matrix()], or write them all with [sw_write()].
#' @param run from [sw_run()]
#' @param what "h" (depth), "eta" (surface elevation, NA where dry),
#'   "speed" (m/s), "hu", "hv"
#' @param dry depth below which a cell counts as dry for "eta" and "speed"
#' @return list of grids, named like `h_t60s`
#' @export
sw_frames <- function(run, what = c("h", "eta", "speed", "hu", "hv"), dry = 1e-3) {
  what <- match.arg(what)
  s <- run$state
  ext <- .run_extent(s)
  out <- lapply(seq_along(run$h), function(k) {
    h <- run$h[[k]]
    m <- switch(what,
      h = h,
      eta = ifelse(h > dry, h + s$z, NA_real_),
      speed = ifelse(h > dry, sqrt(run$hu[[k]]^2 + run$hv[[k]]^2) / h, NA_real_),
      hu = run$hu[[k]], hv = run$hv[[k]])
    .gis(m, ext, s$crs, time = run$diag$time[k])
  })
  names(out) <- sprintf("%s_t%gs", what, run$diag$time)
  out
}

#' Write the frames of a run to a GeoTIFF
#'
#' A multi-band GeoTIFF with one band per output time, the extent and CRS of
#' the state, band descriptions like `h_t60s` and a `time` metadata item
#' (seconds) on each band. Written with gdalraster.
#' @inheritParams sw_frames
#' @param file output file name
#' @param nodata value written for missing cells (dry cells of "eta" and
#'   "speed")
#' @param options GDAL creation options
#' @return `file`, invisibly
#' @export
sw_write <- function(run, file, what = c("h", "eta", "speed", "hu", "hv"), dry = 1e-3,
                     nodata = -9999, options = c("COMPRESS=DEFLATE", "TILED=YES")) {
  .need_gdalraster()
  what <- match.arg(what)
  fr <- sw_frames(run, what, dry)
  s <- run$state
  ext <- .run_extent(s)
  nx <- ncol(s$z); ny <- nrow(s$z)
  ds <- gdalraster::create("GTiff", file, nx, ny, length(fr), "Float64",
                           options = options, return_obj = TRUE)
  on.exit(ds$close())
  ds$setGeoTransform(c(ext[1], (ext[2] - ext[1]) / nx, 0, ext[4], 0, -(ext[4] - ext[3]) / ny))
  wkt <- .as_wkt(s$crs)
  if (nzchar(wkt)) ds$setProjection(wkt)
  for (k in seq_along(fr)) {
    v <- as.vector(fr[[k]])
    v[is.na(v)] <- nodata
    ds$setNoDataValue(k, nodata)
    ds$write(k, 0, 0, nx, ny, v)
    ds$setDescription(k, names(fr)[k])
    ds$setMetadataItem(k, "time", format(attr(fr[[k]], "time"), digits = 15), "")
  }
  invisible(file)
}
