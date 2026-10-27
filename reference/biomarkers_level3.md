# Level 3 (hard) biomarker scenario

6 biomarkers drawn uniformly from the pool (any intensity), small
prevalence gap (~60% vs ~40%), 1.2x intensity multiplier in cases.
Stress-tests the pipeline including dim biomarkers; AUC around 0.70 is a
passing result at this difficulty.

## Usage

``` r
biomarkers_level3(peak_params, H = NULL, W = NULL, seed = NULL)
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
