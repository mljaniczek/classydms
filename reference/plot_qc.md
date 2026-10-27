# Visualize per-sample QC features

Convenience plot showing the cohort distribution of each QC feature,
with flagged samples highlighted. Useful for inspecting whether the
default thresholds are reasonable for your cohort before committing to
exclusion.

## Usage

``` r
plot_qc(
  qc_flagged,
  features = c("total_signal", "mean_positive", "n_peaks", "prop_nonzero",
    "max_intensity", "dynamic_range")
)
```

## Arguments

- qc_flagged:

  Output of
  [`flag_junk_samples()`](https://mljaniczek.github.io/classydms/reference/flag_junk_samples.md).

- features:

  Character vector of features to plot (default: the IQR-checked set).

## Value

A ggplot object (one panel per feature).
