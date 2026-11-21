# Continue from the last frame of a run

Returns a state holding the final depth and momentum of \`run\`, with
any of the forcing fields replaced (for example \`rain = 0\` to stop a
storm). Output times of the next \[sw_run()\] count from zero again.

## Usage

``` r
sw_continue(run, ...)
```

## Arguments

- run:

  from \[sw_run()\]

- ...:

  fields of \[sw_state()\] to replace: manning, rain, f
