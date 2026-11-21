# Set up a shallow water state on a grid

Matrices are oriented like a raster: row 1 is the top (north) edge,
column 1 is the west edge. \`hv\` is momentum towards increasing row
number (southward on a north-up raster).

## Usage

``` r
sw_state(
  z,
  h = 0,
  eta = NULL,
  hu = 0,
  hv = 0,
  dx = 1,
  dy = dx,
  manning = 0,
  rain = 0,
  f = 0,
  g = 9.81,
  extent = NULL,
  crs = NULL
)
```

## Arguments

- z:

  bed elevation matrix (m)

- h:

  initial depth (m), matrix or scalar; or give \`eta\`

- eta:

  optional initial water surface elevation; depth is \`pmax(eta - z,
  0)\`

- hu, hv:

  initial momentum (m^2/s), matrix or scalar

- dx, dy:

  cell size (m)

- manning:

  Manning n (s m^-1/3), matrix or scalar, 0 for frictionless

- rain:

  rain rate (m/s), matrix or scalar

- f:

  Coriolis parameter (1/s), matrix or scalar; 0 for no rotation.
  Positive in the northern hemisphere, negative in the southern; see
  \[sw_coriolis()\]

- g:

  gravity (m/s^2)

- extent, crs:

  optional georeferencing carried through to output
