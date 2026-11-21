# Build a state from a terra DEM

Reads the bed from a single-layer SpatRaster, keeping its extent and CRS
so outputs can be written back as georeferenced rasters. Missing bed
cells are filled with the maximum elevation plus 100 m (a high wall).

## Usage

``` r
sw_from_rast(dem, ...)
```

## Arguments

- dem:

  SpatRaster in metres: projected CRS, or none for a generated bed
  (\[sw_generate()\]); a \`z\` layer is used if there are several
  (\[sw_extract()\])

- ...:

  passed to \[sw_state()\] (h, eta, manning, rain, ...)
