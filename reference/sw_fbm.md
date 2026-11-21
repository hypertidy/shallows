# Fractal noise

Octaves of value noise (fBm), about -1 to 1. Matches playproj
\`src/core/terrain.js\` \`fbm()\`.

## Usage

``` r
sw_fbm(X, Y, scale = 200, octaves = 5, gain = 0.5, seed = 1)
```

## Arguments

- X, Y:

  metres

- scale:

  largest wavelength (m)

- octaves:

  number of octaves

- gain:

  amplitude ratio between octaves

- seed:

  integer seed
