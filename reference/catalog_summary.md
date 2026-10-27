# Diagnostic summary of a peak catalog

Returns a list of scalar and small-data-frame summary stats and, when
`print = TRUE`, prints a human-readable summary. Intended as the first
thing you look at after
[`build_peak_catalog()`](https://mljaniczek.github.io/classydms/reference/build_peak_catalog.md).

## Usage

``` r
catalog_summary(
  catalog,
  prevalence_thresholds = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99),
  print = TRUE
)
```

## Arguments

- catalog:

  A `peak_catalog` object.

- prevalence_thresholds:

  Numeric vector of prevalence values at which to report the number of
  compounds meeting or exceeding that threshold. Default
  `c(0.10, 0.25, 0.50, 0.75, 0.90, 0.99)`.

- print:

  If `TRUE` (default), print a formatted summary in addition to
  returning the underlying data.

## Value

List with `n_pool_peaks`, `n_samples`, `n_compounds`,
`n_small_clusters`, `n_singletons`, `total_observations`,
`fraction_singleton`, `fraction_below_support`, and a data frame
`compounds_at_threshold` with columns
`prevalence_threshold, n_compounds`.
