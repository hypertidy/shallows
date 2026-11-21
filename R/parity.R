## Parity cases shared with the JavaScript port of this package
## (src/core/shallows.js and src/core/terrain.js in hypertidy/playproj).
##
## .parity_spec() lists small problems as plain data. .parity_compute() runs
## them here and attaches the results. inst/parity/make-parity.R writes the
## outcome to inst/parity/parity.json; tests/testthat/test-parity.R checks
## that the package still reproduces that file. The JavaScript tests read a
## copy of the same file and must reproduce it too.
##
## Conventions in the file: grids are row-major with row 0 at the north
## edge; `cells` are 0-based row-major indices (row * nx + col); `bc` is
## west, east, north, south.

.parity_spec <- function() {
  island <- list(shape = "island", nx = 40, ny = 30, dx = 10,
                 noise = list(amp = 3, seed = 7))
  list(
    list(name = "generate_island", kind = "generate", bed = island,
         cells = c(0, 19, 579, 600, 1199)),
    list(name = "generate_reservoir_breach", kind = "generate",
         bed = list(shape = "reservoir", nx = 60, ny = 20, dx = 5,
                    params = list(breach = 20)),
         cells = c(0, 29, 330, 629, 1199)),
    list(name = "fbm_points", kind = "fbm",
         x = c(123.4, 0, -812.25, 5000.5), y = c(-56.7, 0, 77.125, -4321),
         noise = list(scale = 200, octaves = 5, gain = 0.5, seed = 3)),
    list(name = "fill_sea_island", kind = "fill_sea",
         bed = list(shape = "island", nx = 100, ny = 80, dx = 10,
                    noise = list(amp = 3, seed = 7)),
         level = 0),
    list(name = "fill_lake_reservoir", kind = "fill_lake",
         bed = list(shape = "reservoir", nx = 60, ny = 20, dx = 5),
         i0 = 9, j0 = 4, level = 14),
    ## wetting and drying over rough ground, friction and rain, order 2
    list(name = "run_island_slosh", kind = "run", bed = island,
         eta = 20, column = list(i0 = 2, i1 = 6, j0 = 3, j1 = 8, add = 4),
         manning = 0.03, rain = 2e-5, f = 0,
         bc = c("wall", "wall", "wall", "wall"), order = 2, t = 20,
         cells = c(0, 19, 125, 579, 600, 1199)),
    ## open and wall edges with Earth rotation, order 1
    list(name = "run_shelf_rotating", kind = "run",
         bed = list(shape = "shelf", nx = 30, ny = 24, dx = 2000,
                    noise = list(amp = 2, scale = 9000, seed = 11)),
         eta = 0, column = list(i0 = 10, i1 = 13, j0 = 4, j1 = 7, add = 1.5),
         manning = 0.02, rain = 0, f = -1.0e-4,
         bc = c("open", "wall", "wall", "open"), order = 1, t = 1800,
         cells = c(0, 135, 250, 305, 338, 400, 719))
  )
}

.parity_bed <- function(b) {
  args <- c(list(shape = b$shape, nx = b$nx, ny = b$ny, dx = b$dx,
                 noise = b$noise, as = "matrix"), b$params)
  do.call(sw_generate, args)
}

## values of matrix m at 0-based row-major cells
.parity_at <- function(m, cells) as.vector(t(m))[cells + 1L]

.parity_compute <- function(spec = .parity_spec()) {
  lapply(spec, function(s) {
    s$expect <- switch(s$kind,
      generate = {
        z <- .parity_bed(s$bed)
        list(z = .parity_at(z, s$cells), sum = sum(z))
      },
      fbm = list(value = do.call(sw_fbm, c(list(s$x, s$y), s$noise))),
      fill_sea = {
        h <- sw_fill_sea(.parity_bed(s$bed), s$level)
        list(wet = sum(h > 0), sum = sum(h))
      },
      fill_lake = {
        h <- sw_fill_lake(.parity_bed(s$bed), s$i0 + 1L, s$j0 + 1L, s$level)
        list(wet = sum(h > 0), sum = sum(h))
      },
      run = {
        z <- .parity_bed(s$bed)
        h <- pmax(s$eta - z, 0)
        k <- s$column
        rows <- (k$i0:k$i1) + 1L; cols <- (k$j0:k$j1) + 1L
        h[rows, cols] <- h[rows, cols] + k$add
        st <- sw_state(z, h = h, dx = s$bed$dx, manning = s$manning,
                       rain = s$rain, f = s$f)
        run <- sw_run(st, times = s$t, order = s$order, bc = s$bc)
        n <- length(run$h)
        list(h = .parity_at(run$h[[n]], s$cells),
             hu = .parity_at(run$hu[[n]], s$cells),
             hv = .parity_at(run$hv[[n]], s$cells),
             volume = sw_volume(run$h[[n]], st$dx, st$dy),
             wet = sum(run$h[[n]] > 0))
      },
      stop("unknown parity case kind: ", s$kind))
    s
  })
}

.parity_write <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("needs the jsonlite package")
  out <- list(
    package = "shallows",
    version = as.character(utils::packageVersion("shallows")),
    conventions = paste("row-major grids, row 0 north; cells are 0-based",
                        "row * nx + col; bc is west, east, north, south;",
                        "run cases start from depth max(eta - z, 0) plus",
                        "'add' metres on rows i0..i1 and columns j0..j1"),
    cases = .parity_compute())
  writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = I(17), pretty = TRUE), path)
  invisible(path)
}
