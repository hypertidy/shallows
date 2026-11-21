# Frames of a run as a georeferenced SpatRaster

Frames of a run as a georeferenced SpatRaster

## Usage

``` r
sw_rast(run, what = c("h", "eta", "speed", "hu", "hv"), dry = 0.001)
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

SpatRaster with one layer per output time, layer \`time\` set
