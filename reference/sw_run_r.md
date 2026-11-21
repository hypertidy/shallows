# Plain R reference solver (first order)

The same first-order scheme as \[sw_run()\] with \`order = 1\`, written
as vectorised R over whole matrices. It exists to show the method in
readable form and to cross-check the C++ kernel; it is roughly 10 to 50
times slower.

## Usage

``` r
sw_run_r(state, times, cfl = 0.45, bc = "wall", hdry = 1e-05, dt_max = 1)
```

## Arguments

- state:

  from \[sw_state()\]

- times:

  output times (s), increasing; the run ends at the last one

- cfl:

  Courant number; 0.5 is safe for this unsplit scheme

- bc:

  boundary per edge, named \`west\`, \`east\`, \`north\`, \`south\`,
  each "wall", "open" or "periodic"; a single value applies to all.
  Periodic edges come in pairs (west with east, north with south)

- hdry:

  depth below which a cell's velocity is taken as zero

- dt_max:

  largest time step (s); the CFL limit usually decides. A completely dry
  grid steps at 1 s (rain on dry ground)
