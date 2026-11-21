# Build a state from a DEM

Reads the bed from a raster, keeping its extent and CRS so outputs can
be written back as georeferenced rasters. Missing bed cells are filled
with the maximum elevation plus 100 m (a high wall).

## Usage

``` r
sw_from_rast(dem, ...)
```

## Arguments

- dem:

  the bed in metres, in a projected CRS (or none, for a generated bed):
  a file or GDAL source name, a \`gdalraster::GDALRaster\`, or a grid as
  from \[sw_generate()\], \[sw_extract()\] or \`gdalraster::read_ds()\`.
  Only the first band is used

- ...:

  passed to \[sw_state()\] (h, eta, manning, rain, ...)
