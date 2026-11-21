# Package index

## Set up and run

Build a state, run the solver, continue a run.

- [`sw_state()`](https://hypertidy.github.io/shallows/reference/sw_state.md)
  : Set up a shallow water state on a grid
- [`sw_run()`](https://hypertidy.github.io/shallows/reference/sw_run.md)
  : Run the shallow water solver
- [`sw_continue()`](https://hypertidy.github.io/shallows/reference/sw_continue.md)
  : Continue from the last frame of a run
- [`sw_volume()`](https://hypertidy.github.io/shallows/reference/sw_volume.md)
  : Water volume of a depth matrix
- [`sw_coriolis()`](https://hypertidy.github.io/shallows/reference/sw_coriolis.md)
  : Coriolis parameter

## Rasters in and out

Beds from terra SpatRasters, results back as georeferenced rasters.

- [`sw_from_rast()`](https://hypertidy.github.io/shallows/reference/sw_from_rast.md)
  : Build a state from a terra DEM
- [`sw_rast()`](https://hypertidy.github.io/shallows/reference/sw_rast.md)
  : Frames of a run as a georeferenced SpatRaster

## Terrain

Generate synthetic beds, extract real ones, fill lakes and seas.

- [`sw_generate()`](https://hypertidy.github.io/shallows/reference/sw_generate.md)
  : Generate a synthetic bed
- [`sw_shapes()`](https://hypertidy.github.io/shallows/reference/sw_shapes.md)
  : Default parameters of the generator shapes
- [`sw_fbm()`](https://hypertidy.github.io/shallows/reference/sw_fbm.md)
  : Fractal noise
- [`sw_fill_lake()`](https://hypertidy.github.io/shallows/reference/sw_fill_lake.md)
  : Fill a lake
- [`sw_fill_sea()`](https://hypertidy.github.io/shallows/reference/sw_fill_sea.md)
  : Fill the sea
- [`sw_extract()`](https://hypertidy.github.io/shallows/reference/sw_extract.md)
  : Extract a bed for any place
- [`sw_sources()`](https://hypertidy.github.io/shallows/reference/sw_sources.md)
  : Elevation sources known to \[sw_extract()\]
- [`sw_local_crs()`](https://hypertidy.github.io/shallows/reference/sw_local_crs.md)
  : A local metric CRS for a place

## Benchmarks and reference

Analytic solutions used for validation, and the plain R scheme.

- [`dambreak_ritter()`](https://hypertidy.github.io/shallows/reference/dambreak_ritter.md)
  : Ritter dam break onto a dry bed (analytic)
- [`dambreak_stoker()`](https://hypertidy.github.io/shallows/reference/dambreak_stoker.md)
  : Stoker dam break onto a wet bed (analytic)
- [`thacker_bowl()`](https://hypertidy.github.io/shallows/reference/thacker_bowl.md)
  : Thacker planar oscillation in a paraboloid bowl (analytic)
- [`sw_run_r()`](https://hypertidy.github.io/shallows/reference/sw_run_r.md)
  : Plain R reference solver (first order)
