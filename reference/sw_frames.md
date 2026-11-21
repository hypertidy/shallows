# Frames of a run as georeferenced grids

One grid per output time, each a row-major vector with a \`gis\`
attribute (extent and CRS of the state) and a \`time\` attribute, as
read by \`gdalraster::read_ds()\`; plot one with \`ximage::ximage()\`,
get a matrix with \[sw_matrix()\], or write them all with
\[sw_write()\].

## Usage

``` r
sw_frames(run, what = c("h", "eta", "speed", "hu", "hv"), dry = 0.001)
```

## Arguments

- run:

  from \[sw_run()\]

- what:

  "h" (depth), "eta" (surface elevation, NA where dry), "speed" (m/s),
  "hu", "hv"

- dry:

  depth below which a cell counts as dry for "eta" and "speed"

## Value

list of grids, named like \`h_t60s\`
