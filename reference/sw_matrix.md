# Matrix from a grid

A grid as returned by \[sw_generate()\], \[sw_extract()\],
\[sw_frames()\] or \`gdalraster::read_ds()\` is a row-major vector with
a \`gis\` attribute; this gives one band of it as a matrix oriented like
the raster (row 1 at the top), which is the orientation \[sw_state()\]
and \`ximage::ximage()\` use. Anything else is passed to
\[as.matrix()\].

## Usage

``` r
sw_matrix(x, band = 1L)
```

## Arguments

- x:

  grid, or anything \[as.matrix()\] accepts

- band:

  band number

## Value

matrix
