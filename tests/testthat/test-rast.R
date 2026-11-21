test_that("grids round trip through states, frames and GeoTIFF", {
  skip_if_not_installed("gdalraster")
  z <- outer(1:20, 1:30, function(i, j) 0.1 * j)
  ext <- c(500000, 500060, 5200000, 5200040)
  dem <- shallows:::.gis(z, ext, gdalraster::srs_to_wkt("EPSG:28355"))
  expect_equal(sw_matrix(dem), z)
  st <- sw_from_rast(dem, eta = 1.5)
  expect_equal(st$dx, 2)
  expect_equal(st$extent, ext)
  r <- sw_run(st, c(1, 2))
  fr <- sw_frames(r, "h")
  expect_named(fr, c("h_t1s", "h_t2s"))
  expect_equal(sw_matrix(fr[[1]]), r$h[[1]])
  expect_equal(attr(fr[[2]], "gis")$bbox, ext[c(1, 3, 2, 4)])
  expect_equal(attr(fr[[2]], "time"), 2)

  f <- tempfile(fileext = ".tif")
  on.exit(unlink(f))
  sw_write(r, f, "eta")
  ds <- methods::new(gdalraster::GDALRaster, f)
  on.exit(ds$close(), add = TRUE)
  expect_equal(ds$getRasterCount(), 2)
  expect_equal(ds$bbox(), ext[c(1, 3, 2, 4)])
  expect_true(gdalraster::srs_is_same(ds$getProjection(), "EPSG:28355"))
  expect_equal(ds$getDescription(2), "eta_t2s")
  expect_equal(ds$getMetadataItem(2, "time", ""), "2")
  back <- gdalraster::read_ds(ds, bands = 1)
  expect_equal(sw_matrix(back), sw_matrix(sw_frames(r, "eta")[[1]]))

  # and straight back in from the file
  st2 <- sw_from_rast(f)
  expect_equal(st2$extent, ext)
})

test_that("generated grids carry a local extent and plot with ximage", {
  g <- sw_generate("bowl", nx = 40, ny = 20, dx = 10)
  expect_equal(attr(g, "gis")$bbox, c(-200, -100, 200, 100))
  expect_equal(sw_matrix(g), sw_generate("bowl", nx = 40, ny = 20, dx = 10, as = "matrix"),
               ignore_attr = TRUE)
  st <- sw_from_rast(g, eta = 0)
  expect_equal(c(st$dx, st$dy), c(10, 10))
  skip_if_not_installed("ximage")
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({grDevices::dev.off(); unlink(f)})
  expect_equal(ximage::ximage(g)$extent, c(-200, 200, -100, 100))
})

test_that("Coriolis per cell from a state's CRS", {
  skip_if_not_installed("gdalraster")
  g <- shallows:::.gis(matrix(0, 3, 4), c(500000, 540000, 5200000, 5230000), "EPSG:28355")
  f <- sw_coriolis(sw_from_rast(g))
  expect_equal(dim(f), c(3, 4))
  expect_true(all(f < 0))                 # southern hemisphere
  expect_true(all(diff(f[, 1]) < 0))      # row 1 is north: more negative going down the rows
  expect_equal(f, sw_coriolis(g))
})
