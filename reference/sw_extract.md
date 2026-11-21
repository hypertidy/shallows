# Extract a bed for any place

Builds a grid of \`res\` metre cells, \`size\` metres across, centred on
\`lon\`, \`lat\` in a local metric CRS, and fills it from the elevation
sources in order (finest first). Each source is cropped in its own CRS
and warped onto the grid (cell averages when the source is finer than
the grid, bilinear when it is coarser). Later sources fill only cells
the earlier ones left empty. Copernicus GLO-30 has no sea floor (sea
cells are near 0 m), so where GEBCO says sea and GLO-30 is within
\`sea_tol\` of zero, GEBCO's depth is used.

## Usage

``` r
sw_extract(
  lon,
  lat,
  size = 2000,
  res = 10,
  sources = c("tas2m", "glo30", "gebco"),
  crs = sw_local_crs(lon, lat),
  sea_tol = 0.5
)
```

## Arguments

- lon, lat:

  centre (degrees)

- size:

  extent (m), one number for a square or c(width, height)

- res:

  cell size (m)

- sources:

  names from \[sw_sources()\], in the order to use them

- crs:

  target CRS; default \[sw_local_crs()\]

- sea_tol:

  GLO-30 values within this of 0 m count as sea where GEBCO is below 0

## Value

SpatRaster with layers \`z\` (m) and \`source\` (index into \`sources\`,
0 where nothing covered the cell)

## Details

The sources use different vertical datums (see \[sw_sources()\]); they
agree to about a metre around Tasmania, which matters for tides but not
for rain runoff.
