# Choose a common CV window from pooled signal occupancy

Pools per-column signal occupancy (fraction of RT rows above each
sample's robust dust threshold, as in
[`trim_bounds_from_occupancy()`](https://mljaniczek.github.io/classydms/reference/trim_bounds_from_occupancy.md))
across all supplied samples – typically every cohort, excluding any
holdout – on a regular CV axis, and returns the contiguous window of `W`
columns holding the most occupancy. The choice uses no labels.

## Usage

``` r
choose_cv_window(samples, W, dcv = 0.55, cv_anchor = NULL)
```

## Arguments

- samples:

  List of samples with `$cv` and `$Z` (any cohorts; CV sampling may
  differ, each is interpolated onto the common axis).

- W:

  Window width in columns (the model input width).

- dcv:

  CV step (volts) of the common axis.

- cv_anchor:

  A voltage the common axis passes through (keeps columns on the
  instrument's native positions). Default: the most common first CV
  value among the samples.

## Value

A list with `cv` (the `W` window voltages; pass `cv[1]` as `cv0` to
[`make_dms_grid()`](https://mljaniczek.github.io/classydms/reference/make_dms_grid.md)),
`occupancy` (pooled occupancy on the full common axis), `axis` (that
axis), and `captured` (fraction of pooled occupancy inside the window).
