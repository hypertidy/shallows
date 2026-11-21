# Fill a lake

Depth of water standing at \`level\` in the cells connected
(4-neighbour) to cell (\`i0\`, \`j0\`) that lie below it: a reservoir,
crater or hollow full to a level, on any bed. Matches \`fillLake()\` in
playproj \`src/core/terrain.js\`.

## Usage

``` r
sw_fill_lake(z, i0, j0, level)
```

## Arguments

- z:

  bed matrix or grid (see \[sw_matrix()\])

- i0, j0:

  row and column of a cell in the lake (1-based)

- level:

  water level (m)

## Value

depth matrix
