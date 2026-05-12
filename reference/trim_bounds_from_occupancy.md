# Compute occupancy-based trim bounds for one sample

For each RT row and CV column, computes the fraction of pixels above the
dust threshold (the "occupancy"), smooths these curves, and returns the
bounding box of rows/columns whose smoothed occupancy exceeds a
threshold.

## Usage

``` r
trim_bounds_from_occupancy(
  Z,
  eps = NULL,
  thr_rt = 0.005,
  thr_cv = 0.01,
  smooth_k_rt = 31,
  smooth_k_cv = 9,
  min_keep_rt = 100,
  min_keep_cv = 20
)
```

## Arguments

- Z:

  Intensity matrix (post-dust threshold).

- eps:

  Dust threshold (if NULL, computed per-sample as 1% of Q_0.95).

- thr_rt, thr_cv:

  Occupancy thresholds for retaining rows / columns.

- smooth_k_rt, smooth_k_cv:

  Rolling-mean window for occupancy curves.

- min_keep_rt, min_keep_cv:

  Minimum number of rows/columns to keep.

## Value

List with `rt_start`, `rt_end`, `cv_start`, `cv_end`, etc.
