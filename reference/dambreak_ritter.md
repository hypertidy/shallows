# Ritter dam break onto a dry bed (analytic)

Flat frictionless bed, depth \`h0\` for \`x \< x0\`, dry for \`x \>
x0\`.

## Usage

``` r
dambreak_ritter(x, t, h0, x0 = 0, g = 9.81)
```

## Arguments

- x:

  positions (m)

- t:

  time (s), \> 0

- h0:

  upstream depth (m)

- x0:

  dam position

- g:

  gravity

## Value

data frame of x, h, u
