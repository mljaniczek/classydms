# Level 2 (realistic) biomarker scenario

3 biomarkers drawn from the *upper half* of the pool by intensity
(`min_intensity_quantile = 0.50`), moderate prevalence gap (~80% vs
~20%), 1.5x intensity multiplier in cases. Case/control AUC should be at
least ~0.85 if the pipeline works.

## Usage

``` r
biomarkers_level2(peak_params, H = NULL, W = NULL, seed = NULL)
```

## Arguments

- peak_params:

  Output of
  [`estimate_peak_params()`](https://mljaniczek.github.io/classydms/reference/estimate_peak_params.md).

- H, W:

  Optional image dimensions in pixels. When provided, forwarded to
  [`holdout_biomarkers()`](https://mljaniczek.github.io/classydms/reference/holdout_biomarkers.md)
  so the protect radius scales safely for the target image size.
  Strongly recommended — the default fraction-based protect radius
  under-protects on smaller images than ~1400 x 40 and can cause
  background peaks to leak into biomarker regions after location jitter.

- seed:

  Optional RNG seed for reproducible biomarker selection.
