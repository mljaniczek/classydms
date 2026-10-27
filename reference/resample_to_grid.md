# Resample one sample onto a dms_grid (RT and CV)

Applies
[`resample_to_time_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_time_grid.md)
with `grid$time`, then
[`resample_to_cv_grid()`](https://mljaniczek.github.io/classydms/reference/resample_to_cv_grid.md)
with `grid$cv`. The result is exactly `grid$H x grid$W`, so every sample
shares one physical pixel grid and needs no padding.

## Usage

``` r
resample_to_grid(sample, grid, method = c("linear", "constant"))
```

## Arguments

- sample:

  A sample list with `$time`, `$cv`, `$Z`.

- grid:

  A
  [`dms_grid()`](https://mljaniczek.github.io/classydms/reference/dms_grid.md).

- method:

  Interpolation method for both axes.

## Value

The resampled sample (`$time = grid$time`, `$cv = grid$cv`).
