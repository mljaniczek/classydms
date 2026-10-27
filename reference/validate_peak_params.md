# Validate the internal consistency of a peak_params object

Checks the invariants every function that consumes `peak_params` relies
on:

- all per-peak vectors (every field ending in `_raw`, e.g. `rt_loc_raw`,
  `rt_seconds_raw`, `sample_idx_raw`) have one common length;

- `rt_seconds_raw` and `cv_volts_raw` are present together or not at
  all;

- `sample_idx_raw`, if present, is integer-valued with no `NA` and lies
  in `1..n_samples`;

- per-sample fields (`sample_names`, `sample_cohort`, `n_peaks$values`,
  the per-sample vectors in `alignment`) have length `n_samples`, and
  `time_ranges` / `cv_ranges` have `n_samples` rows;

- `rt_loc$values` / `cv_loc$values`, if present, match the per-peak
  length.

Called at the top of
[`build_peak_catalog()`](https://mljaniczek.github.io/classydms/reference/build_peak_catalog.md)
and at the end of every function that returns a transformed
`peak_params`. Useful on saved objects too: a physical-mode
`peak_params` written by
[`subset_peak_params()`](https://mljaniczek.github.io/classydms/reference/subset_peak_params.md)
before this check existed fails here (its `rt_seconds_raw` /
`cv_volts_raw` were left unsubset).

## Usage

``` r
validate_peak_params(peak_params, context = "peak_params")
```

## Arguments

- peak_params:

  A `peak_params` list.

- context:

  Label used in error messages.

## Value

`peak_params`, invisibly. Errors on the first violation.
