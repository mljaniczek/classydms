# Compute a common regular time grid across a cohort of samples

Chooses a regular RT grid suitable for resampling every sample onto
shared row indices. Row `i` in the resampled Z always corresponds to
physical time `target_time[i]` regardless of which sample the row came
from — the prerequisite for pixel-consistent cross-sample encoder input.

## Usage

``` r
compute_common_time_grid(samples, step = 1, t_max = NULL)
```

## Arguments

- samples:

  List of samples (each with a `$time` numeric vector).

- step:

  Grid step in seconds. See details.

- t_max:

  Maximum RT for the grid (seconds). See details.

## Value

Numeric vector of RT positions in seconds, evenly spaced by `step`,
starting at 0, ending at `>= max sample time`.

## Details

Grid choice:

- `step`: distance in seconds between consecutive rows. Default `1.0` —
  physically meaningful (row index = seconds), adequate for typical
  GC-DMS peak FWHM (~2-5 s), and roughly matches current native RT
  dimensions (~1500 rows). Common alternatives: `0.5` (double
  resolution), `median(diff(time))` from the cohort, or `"median"` /
  `"mean"` which computes the median / mean per-row time delta across
  all samples.

- `t_max`: maximum RT (seconds) on the returned grid. Default `NULL`
  uses `max(unlist(lapply(samples, `\[\[`, "time")))` — the longest
  sample's endpoint, so no sample is truncated. Rounded up to the
  nearest integer multiple of `step`.
