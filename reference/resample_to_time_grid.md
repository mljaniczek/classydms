# Resample one sample onto a regular time grid via linear interpolation

Interpolates each CV column of the sample's `$Z` matrix from the
sample's own `$time` vector onto a target regular grid. Values past the
sample's own maximum RT are filled with 0 (no signal). Returns a new
sample with `$Z` on the target grid and `$time` set to the target grid.

## Usage

``` r
resample_to_time_grid(sample, target_time, method = c("linear", "constant"))
```

## Arguments

- sample:

  A sample list with `$time` (numeric vector) and `$Z` (matrix with
  `nrow(Z) == length(time)`).

- target_time:

  Numeric vector of RT positions (seconds) — the grid to interpolate
  onto. Typically from
  [`compute_common_time_grid()`](https://mljaniczek.github.io/classydms/reference/compute_common_time_grid.md).

- method:

  Interpolation method passed to
  [`stats::approx()`](https://rdrr.io/r/stats/approxfun.html). Default
  `"linear"`. Alternatives: `"constant"` (nearest-neighbor step
  function).

## Value

A sample list with updated `$Z` (dims `length(target_time) x ncol(Z)`)
and `$time` set to `target_time`.

## Details

Purpose: after resampling, row index `i` in every sample corresponds to
the same physical RT (`target_time[i]`), removing cross-sample
misalignment caused by variable per-row time deltas and different RT
endpoints. Downstream encoders / classifiers can then treat pixel
positions as physically meaningful and cross-sample-comparable.

`$cv` is unchanged (CV grid is already consistent across samples on the
same instrument). Other list fields (`$path`, `$id`, etc.) are
preserved.
