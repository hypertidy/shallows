# Coriolis parameter

\`2 \* omega \* sin(latitude)\`: positive in the northern hemisphere and
negative in the southern. Give a state or a grid with a CRS to get one
value per cell (cell centres are transformed to longitude and latitude
with gdalraster), or a numeric vector of latitudes in degrees.

## Usage

``` r
sw_coriolis(x, omega = 7.2921159e-05)
```

## Arguments

- x:

  a state from \[sw_state()\] or \[sw_from_rast()\] with extent and CRS,
  a grid with a CRS (see \[sw_matrix()\]), or latitudes in degrees

- omega:

  Earth's rotation rate (rad/s)

## Value

matrix (for a state or grid, oriented like the bed) or vector
