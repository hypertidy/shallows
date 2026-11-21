# Generate a synthetic bed

Shapes you can reason about, with optional fractal roughness on top.
Local coordinates are metres from the grid centre, X east and Y north;
row 1 is the north edge, as in a raster.

## Usage

``` r
sw_generate(
  shape = names(.shapes),
  nx = 200,
  ny = nx,
  dx = 5,
  dy = dx,
  ...,
  noise = NULL,
  as = c("grid", "matrix")
)
```

## Arguments

- shape:

  one of the names above

- nx, ny:

  cells across and down

- dx, dy:

  cell size (m)

- ...:

  shape parameters, replacing the defaults (see \`sw_shapes()\`)

- noise:

  NULL, or a list with \`amp\` (m) and optionally \`scale\`,
  \`octaves\`, \`gain\`, \`seed\` for \[sw_fbm()\]

- as:

  "grid" for a grid in local metres with no CRS (a row-major vector with
  a \`gis\` attribute, as from \`gdalraster::read_ds()\`; plot it with
  \`ximage::ximage()\`, give it to \[sw_from_rast()\]) or "matrix" (with
  \`dx\`, \`dy\` attributes)

## Details

\* \`plane\`: \`z0 - slope_x X - slope_y Y\`. \* \`valley\`: drains
east, floor falling \`slope\`, sides rising \`side\` beyond a flat floor
\`floor\` wide. \* \`bowl\`: paraboloid, \`depth\` below rim level at
the centre, rim at \`radius\`. \* \`shelf\`: north-south coast at X = 0,
land west rising \`land_slope\`, shelf east falling \`shelf_slope\` to
\`shelf_depth\`, then \`slope_slope\` to \`abyss\`. \* \`island\`:
gaussian \`height\` above sea level on a floor \`depth\` deep. \*
\`reservoir\`: a valley with a dam at \`dam_x\`, \`dam_height\` above
the floor there, \`dam_thick\` thick, breached over \`breach\` metres (0
= intact).
