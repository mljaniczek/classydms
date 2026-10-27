# Compute per-sample QC features

Calculates a battery of sample-level summary statistics on a list of
preprocessed sample matrices (typically after Savitzky-Golay + AsLS +
dust thresholding + trim, but before normalization or padding). Returns
a tibble with one row per sample and columns covering signal intensity,
sparsity, peak density, and (optional) cohort-relative similarity.

## Usage

``` r
qc_features(samples, sample_names = NULL, compute_cohort_correlation = FALSE)
```

## Arguments

- samples:

  A list of sample objects, each with element `Z` (the preprocessed
  intensity matrix). Typically the output of
  `lapply(samples, trim_sample, bounds = global_b)`.

- sample_names:

  Optional character vector of sample names. If `NULL`, uses
  `names(samples)` or falls back to integer indices.

- compute_cohort_correlation:

  If `TRUE`, also compute each sample's correlation against the cohort
  median spectrum. Slower (O(n_samples \* n_pixels)) but useful for
  detecting morphological outliers that don't trigger univariate flags.
  Default `FALSE`.

## Value

A tibble with one row per sample and columns:

- `sample` — sample identifier

- `rt_dim`, `cv_dim` — post-trim dimensions

- `n_positive`, `prop_nonzero` — count and fraction of pixels above zero

- `total_signal` — sum of all intensities

- `mean_positive` — mean of pixels above zero

- `max_intensity` — peak pixel value

- `dynamic_range` — `max_intensity / median(positive pixels)`

- `bbox_size` — `rt_dim * cv_dim`

- `n_peaks` — detected peaks via
  [`pick_peak_centers()`](https://mljaniczek.github.io/classydms/reference/pick_peak_centers.md)

- `cohort_correlation` — (if requested) correlation with cohort median

## Details

The returned table is consumed by
[`flag_junk_samples()`](https://mljaniczek.github.io/classydms/reference/flag_junk_samples.md),
which applies user-tunable thresholds to mark samples for exclusion.
