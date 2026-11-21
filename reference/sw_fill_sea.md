# Fill the sea

Depth of water standing at \`level\` in cells below it that connect
(4-neighbour) to the edge of the grid; inland hollows stay dry. Matches
\`fillSea()\` in playproj \`src/core/terrain.js\`.

## Usage

``` r
sw_fill_sea(z, level = 0)
```

## Arguments

- z:

  bed matrix

- level:

  sea level (m)
