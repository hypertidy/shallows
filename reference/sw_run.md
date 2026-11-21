# Run the shallow water solver

Run the shallow water solver

## Usage

``` r
sw_run(
  state,
  times,
  order = 2L,
  cfl = 0.45,
  bc = "wall",
  hdry = 1e-05,
  max_steps = 1e+07,
  dt_max = Inf
)
```

## Arguments

- state:

  from \[sw_state()\]

- times:

  output times (s), increasing; the run ends at the last one

- order:

  1 (first order) or 2 (MUSCL + Heun)

- cfl:

  Courant number; 0.5 is safe for this unsplit scheme

- bc:

  boundary per edge, named \`west\`, \`east\`, \`north\`, \`south\`,
  each "wall", "open" or "periodic"; a single value applies to all.
  Periodic edges come in pairs (west with east, north with south)

- hdry:

  depth below which a cell's velocity is taken as zero

- max_steps:

  safety limit on time steps

- dt_max:

  largest time step (s); the CFL limit usually decides. A completely dry
  grid steps at 1 s (rain on dry ground)

## Value

list with \`h\`, \`hu\`, \`hv\` (lists of matrices, one per time),
\`diag\` (data frame: volume, cumulative outflow and rain, budget error)
and the input \`state\`
