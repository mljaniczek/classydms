# Subset a peak_params object to a chosen set of samples

Given a `peak_params` object with `sample_idx_raw` (or after
[`backfill_sample_idx_raw()`](https://mljaniczek.github.io/classydms/reference/backfill_sample_idx_raw.md)),
returns a new `peak_params` containing only the peaks and metadata for
the specified samples. Useful for building per-condition catalogs from a
single pooled `estimate_peak_params` run — split a combined
disease+control `peak_params` into `pp_disease` and `pp_control` without
rerunning the expensive step.

## Usage

``` r
subset_peak_params(peak_params, sample_indices)
```

## Arguments

- peak_params:

  A `peak_params` object with `sample_idx_raw`.

- sample_indices:

  Integer vector of sample indices (in the original peak_params'
  1..n_samples space) to keep. Common pattern:
  `which(y == disease_label)`.

## Value

A `peak_params`-shaped list containing only the specified samples'
peaks. Ready to pass to
[`build_peak_catalog()`](https://mljaniczek.github.io/classydms/reference/build_peak_catalog.md)
and other downstream tools.

## Details

Sample indices in the returned object are remapped to
`1:length(sample_indices)` so downstream tools that expect contiguous
1-based indexing work. Samples keep their original relative order
(`sample_indices` is sorted). Every per-peak vector (every field ending
in `_raw`, including `rt_seconds_raw` / `cv_volts_raw`) is subset with
the same mask, and every per-sample field (`sample_names`,
`sample_cohort`, `time_ranges`, `cv_ranges`, `n_peaks$values`,
per-sample `alignment` vectors) is subset to the kept samples. Fitted
summary statistics (`rt_loc$mean/sd`, `sigma_rt$meanlog/sdlog`, etc.)
are recomputed from the subset; noise stats are kept as-is (pooled
across samples in the original run, which is usually close enough unless
the condition of interest has systematically different background
noise).
