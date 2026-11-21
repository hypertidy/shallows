# Write the frames of a run to a GeoTIFF

A multi-band GeoTIFF with one band per output time, the extent and CRS
of the state, band descriptions like \`h_t60s\` and a \`time\` metadata
item (seconds) on each band. Written with gdalraster.

## Usage

``` r
sw_write(
  run,
  file,
  what = c("h", "eta", "speed", "hu", "hv"),
  dry = 0.001,
  nodata = -9999,
  options = c("COMPRESS=DEFLATE", "TILED=YES")
)
```

## Arguments

- run:

  from \[sw_run()\]

- file:

  output file name

- what:

  "h" (depth), "eta" (surface elevation, NA where dry), "speed" (m/s),
  "hu", "hv"

- dry:

  depth below which a cell counts as dry for "eta" and "speed"

- nodata:

  value written for missing cells (dry cells of "eta" and "speed")

- options:

  GDAL creation options

## Value

\`file\`, invisibly
