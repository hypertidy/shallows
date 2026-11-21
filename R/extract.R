#' Elevation sources known to [sw_extract()]
#'
#' A data frame of the sources, finest first: name, nominal cell size (m),
#' vertical datum, coverage test, and how to open the window. Sources are
#' read in place from cloud-optimised GeoTIFFs through GDAL `/vsicurl/`.
#' @export
sw_sources <- function() {
  data.frame(
    name = c("tas2m", "glo30", "gebco"),
    res_m = c(2, 30, 460),
    datum = c("AHD (Australian Height Datum)", "EGM2008 geoid", "mean sea level"),
    covers = c("Tasmania, land", "global land (no ocean depths)", "global, land and sea floor"),
    url = c("https://data.source.coop/alexgleith/tasmania-dem-2m/Tasmania_Statewide_2m_DEM_14-08-2021.tif",
            "https://copernicus-dem-30m.s3.amazonaws.com/ (1 degree tiles)",
            "https://data.source.coop/ausantarctic/gebco/GEBCO_2026.tif"),
    stringsAsFactors = FALSE)
}

#' A local metric CRS for a place
#'
#' UTM (WGS84) between 80 S and 84 N, polar stereographic beyond.
#' @param lon,lat centre in degrees
#' @export
sw_local_crs <- function(lon, lat) {
  if (lat < -80) return("EPSG:3031")
  if (lat > 84) return("EPSG:3413")
  zone <- floor((lon + 180) / 6) %% 60 + 1
  sprintf("EPSG:%d", (if (lat < 0) 32700 else 32600) + zone)
}

.glo30_tiles <- function(ll_ext) {
  lats <- seq(floor(ll_ext[3]), floor(ll_ext[4] - 1e-9))
  lons <- seq(floor(ll_ext[1]), floor(ll_ext[2] - 1e-9))
  g <- expand.grid(lat = lats, lon = lons)
  nm <- sprintf("Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM",
                ifelse(g$lat < 0, "S", "N"), abs(g$lat), ifelse(g$lon < 0, "W", "E"), abs(g$lon))
  sprintf("/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/%s/%s.tif", nm, nm)
}

# GDAL source names to warp from, or NULL when the source cannot cover the
# window (ll_ext is xmin, xmax, ymin, ymax in degrees)
.open_source <- function(name, ll_ext) {
  if (name == "tas2m") {
    # Tasmania and its islands, roughly
    if (ll_ext[2] < 143.5 || ll_ext[1] > 148.6 || ll_ext[4] < -43.8 || ll_ext[3] > -39.4) return(NULL)
    return(paste0("/vsicurl/", sw_sources()$url[1]))
  }
  if (name == "glo30") {
    tiles <- .glo30_tiles(ll_ext)
    ok <- tiles[vapply(tiles, function(u) isTRUE(gdalraster::vsi_stat(u)), logical(1))]
    if (!length(ok)) return(NULL)  # open ocean: no tiles
    return(ok)
  }
  if (name == "gebco") return(paste0("/vsicurl/", sw_sources()$url[3]))
  stop("unknown source ", name)
}

# Warp source(s) onto the target grid (te = xmin, ymin, xmax, ymax in the
# target CRS); row-major values with NA where the source has no data, or
# NULL if a single source does not reach the window
.warp_source <- function(dsn, te, nx, ny, wkt, res) {
  ds <- methods::new(gdalraster::GDALRaster, dsn[1])
  srs <- ds$getProjection(); r <- ds$res(); bb <- ds$bbox()
  ds$close()
  if (length(dsn) == 1L) {
    # window in the source's CRS, with a margin of a few source cells
    w <- gdalraster::transform_bounds(te, wkt, srs)
    m <- 4 * max(r)
    if (w[1] > bb[3] + m || w[3] < bb[1] - m || w[2] > bb[4] + m || w[4] < bb[2] - m) return(NULL)
  }
  finer <- max(r) * (if (gdalraster::srs_is_geographic(srs)) 111000 else 1) < res
  out <- paste0("/vsimem/", basename(tempfile("shallows_")), ".tif")
  on.exit(gdalraster::vsi_unlink(out))
  # gdalwarp reads the right overview level of a COG, so a coarse grid over
  # a fine source reads only what it needs
  gdalraster::warp(dsn, out, wkt,
                   cl_arg = c("-te", sprintf("%.17g", te), "-ts", nx, ny,
                              "-r", if (finer) "average" else "bilinear",
                              "-ot", "Float64", "-dstnodata", "-99999"),
                   quiet = TRUE)
  o <- methods::new(gdalraster::GDALRaster, out)
  v <- as.vector(gdalraster::read_ds(o))
  o$close()
  v[!is.finite(v) | v == -99999] <- NA
  v
}

#' Extract a bed for any place
#'
#' Builds a grid of `res` metre cells, `size` metres across, centred on
#' `lon`, `lat` in a local metric CRS, and fills it from the elevation
#' sources in order (finest first). Each source is warped onto the grid with
#' gdalraster (cell averages when the source is finer than the grid,
#' bilinear when it is coarser). Later sources fill only cells the earlier
#' ones left empty. Copernicus GLO-30 has no sea floor (sea cells are near
#' 0 m), so where GEBCO says sea and GLO-30 is within `sea_tol` of zero,
#' GEBCO's depth is used.
#'
#' The sources use different vertical datums (see [sw_sources()]); they
#' agree to about a metre around Tasmania, which matters for tides but not
#' for rain runoff.
#'
#' @param lon,lat centre (degrees)
#' @param size extent (m), one number for a square or c(width, height)
#' @param res cell size (m)
#' @param sources names from [sw_sources()], in the order to use them
#' @param crs target CRS; default [sw_local_crs()]
#' @param sea_tol GLO-30 values within this of 0 m count as sea where GEBCO
#'   is below 0
#' @return the bed (m) as a grid: a row-major vector with a `gis` attribute
#'   (extent and CRS, as from `gdalraster::read_ds()`), ready for
#'   [sw_from_rast()] and `ximage::ximage()`. Attribute `source` is an
#'   integer matrix saying which of `sources` filled each cell (0 where
#'   none did); attributes `sources`, `lon` and `lat` record the call
#' @export
sw_extract <- function(lon, lat, size = 2000, res = 10,
                       sources = c("tas2m", "glo30", "gebco"),
                       crs = sw_local_crs(lon, lat), sea_tol = 0.5) {
  .need_gdalraster()
  size <- rep(size, length.out = 2)
  wkt <- gdalraster::srs_to_wkt(crs)
  ctr <- gdalraster::transform_xy(cbind(lon, lat), "EPSG:4326", wkt)
  nx <- max(2L, round(size[1] / res)); ny <- max(2L, round(size[2] / res))
  te <- c(ctr[1] - nx * res / 2, ctr[2] - ny * res / 2, ctr[1] + nx * res / 2, ctr[2] + ny * res / 2)
  ll <- gdalraster::transform_bounds(te, wkt, "EPSG:4326")[c(1, 3, 2, 4)]
  z <- rep(NA_real_, nx * ny); src <- integer(nx * ny)
  gebco <- NULL
  for (k in seq_along(sources)) {
    dsn <- .open_source(sources[k], ll)
    if (is.null(dsn)) next
    rw <- .warp_source(dsn, te, nx, ny, wkt, res)
    if (is.null(rw)) next
    if (sources[k] == "gebco") gebco <- rw
    fill <- is.na(z) & !is.na(rw)
    z[fill] <- rw[fill]
    src[fill] <- k
  }
  # GLO-30 sea cells: take GEBCO's depth
  kg <- match("glo30", sources)
  if (!is.na(kg) && !is.null(gebco)) {
    sea <- src == kg & !is.na(z) & abs(z) <= sea_tol & !is.na(gebco) & gebco < 0
    z[sea] <- gebco[sea]
    src[sea] <- match("gebco", sources)
  }
  .gis(matrix(z, ny, nx, byrow = TRUE), te[c(1, 3, 2, 4)], wkt,
       source = matrix(src, ny, nx, byrow = TRUE), sources = sources, lon = lon, lat = lat)
}
