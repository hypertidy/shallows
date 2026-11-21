# Coriolis parameter

\`2 \* omega \* sin(latitude)\`: positive in the northern hemisphere and
negative in the southern. Give a SpatRaster to get one value per cell
(its cell centres are transformed to longitude and latitude), or a
numeric vector of latitudes in degrees.

## Usage

``` r
sw_coriolis(x, omega = 7.2921159e-05)
```

## Arguments

- x:

  SpatRaster, or latitudes in degrees

- omega:

  Earth's rotation rate (rad/s)

## Value

matrix (for a raster, oriented like \[sw_from_rast()\]) or vector
