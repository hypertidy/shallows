# Thacker planar oscillation in a paraboloid bowl (analytic)

Bed \`z = D0 \* (r^2 / L^2 - 1)\` centred at the origin. The water
surface stays a plane that tilts and circles around the bowl with period
\`2 pi / omega\`, \`omega = sqrt(2 g D0) / L\`, and the shoreline moves:
a standard wetting and drying benchmark (Thacker 1981, J. Fluid Mech.
107, 499-508).

## Usage

``` r
thacker_bowl(x, y, t, D0 = 0.1, L = 1, eta0 = 0.5, g = 9.81)
```

## Arguments

- x, y:

  positions (m), recycled to common length

- t:

  time (s)

- D0:

  depth at the centre at rest (m)

- L:

  radius of the rest shoreline (m)

- eta0:

  displacement of the shoreline circle's centre (m); the defaults are
  the SWASHES benchmark (a = 1, h0 = 0.1, eta = 0.5)

- g:

  gravity

## Value

list with bed \`z\`, depth \`h\`, velocities \`u\`, \`v\`, and
\`period\`
