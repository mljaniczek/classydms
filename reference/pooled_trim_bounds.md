# Pool per-sample trim bounds to a single cohort-wide bounding box

Uses robust quantiles (5th of starts, 95th of ends by default) so a
single outlier sample with unusually narrow or wide content does not
force the cohort-wide bounds to lose data common to nearly all samples.

## Usage

``` r
pooled_trim_bounds(bounds_list, R, C, q_lo = 0.05, q_hi = 0.95)
```

## Arguments

- bounds_list:

  List of bounds (one per sample) from
  [`trim_bounds_from_occupancy()`](https://mljaniczek.github.io/classydms/reference/trim_bounds_from_occupancy.md).

- R, C:

  Maximum row / column counts across samples; the pooled bounds are
  clamped to lie within `[1, R]` and `[1, C]`.

- q_lo, q_hi:

  Quantiles for pooling start / end positions (default 0.05 and 0.95).

## Value

List with `rt_start`, `rt_end`, `cv_start`, `cv_end` giving the
cohort-wide bounding box.
