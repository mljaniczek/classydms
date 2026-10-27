# Plot real vs synthetic per-sample summary statistics

Faceted ggplot showing the distribution of each per-sample summary
statistic (total signal, mean positive intensity, sparsity, peak count)
for real and synthetic samples side-by-side. Use as a visual companion
to
[`synthetic_quality_check()`](https://mljaniczek.github.io/classydms/reference/synthetic_quality_check.md).

## Usage

``` r
plot_synthetic_quality(qc_result)
```

## Arguments

- qc_result:

  Output of
  [`synthetic_quality_check()`](https://mljaniczek.github.io/classydms/reference/synthetic_quality_check.md).

## Value

A ggplot object.
